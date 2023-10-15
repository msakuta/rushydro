use std::{
    cell::Cell,
    collections::{HashMap, VecDeque},
};

use eframe::{
    egui::{
        self,
        plot::{Legend, Line, Plot, PlotPoints},
        Context, Frame, Ui,
    },
    emath::Align2,
    epaint::{pos2, vec2, Color32, FontId, Pos2, Rect, Vec2},
};
use rand::{thread_rng, Rng};

use crate::measure_time;

const SCALE: f32 = 10.;
const NUM_PARTICLES: usize = 200;
const PARTICLE_RADIUS: f32 = 2.;
const PARTICLE_RADIUS2: f32 = PARTICLE_RADIUS * PARTICLE_RADIUS;
const PARTICLE_RENDER_RADIUS: f32 = 3.;
const CELL_SIZE: f32 = PARTICLE_RADIUS;
const RESTITUTION: f32 = 0.5;
const REPULSION_FORCE: f32 = 0.1;
const VISCOSITY: f32 = 0.01;
const SURFACE_TENSION: f32 = 0.05;
const G: f32 = 0.01;
const MOUSE_EFFECT_RADIUS: f32 = 5.;
const MAX_HISTORY: usize = 1000;

struct Particle {
    pos: Cell<Vec2>,
    velo: Cell<Vec2>,
}

#[derive(PartialEq, Eq, Debug)]
enum NeighborMode {
    BruteForce,
    HashMap,
    SortMap,
}

enum NeighborPayload {
    BruteForce,
    HashMap(HashMap<(i32, i32), Vec<usize>>),
    SortMap {
        hash_table: Vec<HashEntry>,
        start_offsets: Vec<usize>,
    },
}

#[derive(Default, Clone, Copy, Debug)]
struct HashEntry {
    particle_idx: usize,
    cell_hash: usize,
}

impl std::fmt::Display for HashEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{{}, {}}}", self.particle_idx, self.cell_hash)
    }
}

struct Params {
    repulsion_force: f32,
    viscosity: f32,
    surface_tension: f32,
}

pub struct RusHydroApp {
    particles: Vec<Particle>,
    rect: Rect,
    num_particles: usize,
    restitution: f32,
    params: Params,
    gravity: f32,
    mouse_down: Option<(Pos2, Vec2)>,
    mouse_effect_radius: f32,
    neighbor_mode: NeighborMode,
    neighbor_payload: NeighborPayload,
    show_grid: bool,
    show_grid_count: bool,
    color_by_speed: bool,
    show_time_plot: bool,
    time_history: VecDeque<f64>,
    hash_time_history: VecDeque<f64>,
    sort_time_history: VecDeque<(f64, f64)>,
}

impl RusHydroApp {
    pub fn new() -> Self {
        let mut rng = thread_rng();
        let rect = Rect::from_center_size(Pos2::ZERO, Vec2::splat(30.));
        let particles = (0..NUM_PARTICLES)
            .map(|_| Particle {
                pos: Cell::new(vec2(
                    rng.gen_range(rect.min.x..rect.max.x),
                    rng.gen_range(rect.min.y..rect.max.y),
                )),
                velo: Cell::new(Vec2::ZERO),
            })
            .collect();
        Self {
            particles,
            rect,
            num_particles: NUM_PARTICLES,
            restitution: RESTITUTION,
            params: Params {
                repulsion_force: REPULSION_FORCE,
                viscosity: VISCOSITY,
                surface_tension: SURFACE_TENSION,
            },
            gravity: G,
            mouse_down: None,
            mouse_effect_radius: MOUSE_EFFECT_RADIUS,
            neighbor_mode: NeighborMode::HashMap,
            neighbor_payload: NeighborPayload::BruteForce,
            show_grid: false,
            show_grid_count: false,
            color_by_speed: true,
            show_time_plot: false,
            time_history: VecDeque::new(),
            hash_time_history: VecDeque::new(),
            sort_time_history: VecDeque::new(),
        }
    }

    fn paint_canvas(&mut self, ui: &mut Ui) {
        Frame::canvas(ui.style()).show(ui, |ui| {
            let (response, painter) =
                ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());

            let to_screen = egui::emath::RectTransform::from_to(
                Rect::from_min_size(Pos2::ZERO, response.rect.size()),
                response.rect,
            );
            let from_screen = to_screen.inverse();

            let canvas_offset_x = response.rect.width() * 0.5;
            let canvas_offset_y = response.rect.height() * 0.5;

            let to_pos2 = |pos: Pos2| {
                to_screen.transform_pos(pos2(
                    canvas_offset_x + pos.x as f32 * SCALE,
                    canvas_offset_y - pos.y as f32 * SCALE,
                ))
            };

            let from_pos2 = |pos: Pos2| {
                let scr_pos = from_screen.transform_pos(pos);
                pos2(
                    (scr_pos.x - canvas_offset_x) / SCALE,
                    -(scr_pos.y - canvas_offset_y) / SCALE,
                )
            };

            match self.neighbor_payload {
                NeighborPayload::HashMap(ref hash_map) => {
                    if self.show_grid {
                        let font_id = FontId::monospace(16.);
                        for (grid_pos, cell) in hash_map {
                            let rect = Rect::from_min_size(
                                pos2(grid_pos.0 as f32 * CELL_SIZE, grid_pos.1 as f32 * CELL_SIZE),
                                Vec2::splat(CELL_SIZE),
                            );
                            let mut scr_rect =
                                Rect::from_min_max(to_pos2(rect.min), to_pos2(rect.max));
                            std::mem::swap(&mut scr_rect.min.y, &mut scr_rect.max.y);
                            painter.rect_stroke(
                                scr_rect,
                                0.,
                                (1., Color32::from_rgb(255, 127, 127)),
                            );
                            if self.show_grid_count {
                                painter.text(
                                    scr_rect.min,
                                    Align2::LEFT_TOP,
                                    format!("{}", cell.len()),
                                    font_id.clone(),
                                    Color32::BLACK,
                                );
                            }
                        }
                    }
                }
                NeighborPayload::SortMap { .. } => {
                    if self.show_grid {
                        for particle in &self.particles {
                            let pos = particle.pos.get();
                            let grid_pos = (
                                pos.x.div_euclid(CELL_SIZE) as i32,
                                pos.y.div_euclid(CELL_SIZE) as i32,
                            );
                            let rect = Rect::from_min_size(
                                pos2(grid_pos.0 as f32 * CELL_SIZE, grid_pos.1 as f32 * CELL_SIZE),
                                Vec2::splat(CELL_SIZE),
                            );
                            let mut scr_rect =
                                Rect::from_min_max(to_pos2(rect.min), to_pos2(rect.max));
                            std::mem::swap(&mut scr_rect.min.y, &mut scr_rect.max.y);
                            painter.rect_stroke(
                                scr_rect,
                                0.,
                                (1., Color32::from_rgb(255, 127, 127)),
                            );
                        }
                    }
                }
                _ => {}
            }

            let mut scr_rect = Rect::from_min_max(to_pos2(self.rect.min), to_pos2(self.rect.max));
            std::mem::swap(&mut scr_rect.min.y, &mut scr_rect.max.y);
            scr_rect = scr_rect.expand(PARTICLE_RENDER_RADIUS);
            painter.rect_stroke(scr_rect, 0., (1., Color32::BLACK));

            if response.is_pointer_button_down_on() {
                if let Some(ptr) = response.interact_pointer_pos() {
                    let pos = from_pos2(ptr);
                    let delta = self
                        .mouse_down
                        .map(|(prev_pos, _)| pos - prev_pos)
                        .unwrap_or(Vec2::ZERO);
                    self.mouse_down = Some((pos, delta));

                    painter.circle(to_pos2(pos), 5., Color32::RED, (1., Color32::BLACK));
                } else {
                    self.mouse_down = None;
                }
            } else {
                self.mouse_down = None;
            }

            for particle in &self.particles {
                let color = if self.color_by_speed {
                    let speed = particle.velo.get().length().min(1.);
                    let speed_u8 = (speed * 255.) as u8;
                    Color32::from_rgb(speed_u8, 0, 255 - speed_u8)
                } else {
                    Color32::BLUE
                };
                painter.circle_filled(
                    to_pos2(particle.pos.get().to_pos2()),
                    PARTICLE_RENDER_RADIUS,
                    color,
                );
            }
        });
    }

    fn reset(&mut self) {
        let mut rng = thread_rng();
        let particles = (0..self.num_particles)
            .map(|_| Particle {
                pos: Cell::new(vec2(
                    rng.gen_range(self.rect.min.x..self.rect.max.x),
                    rng.gen_range(self.rect.min.y..self.rect.max.y),
                )),
                velo: Cell::new(Vec2::ZERO),
            })
            .collect();
        self.particles = particles;
        // Force reinitialization of neighbor caches
        self.neighbor_payload = NeighborPayload::BruteForce;
    }

    fn update_speed(particle_i: &Particle, particle_j: &Particle, params: &Params) {
        let pos_i = particle_i.pos.get();
        let pos_j = particle_j.pos.get();
        let delta = pos_i - pos_j;
        let velo_i = particle_i.velo.get();
        let velo_j = particle_j.velo.get();
        let average_velo = (velo_i + velo_j) * 0.5;
        let dist2 = delta.length_sq();
        if 0. < dist2 && dist2 < PARTICLE_RADIUS2 {
            let dist = dist2.sqrt();
            let repulsion =
                delta / dist * ((1. - dist / PARTICLE_RADIUS).powf(2.) - params.surface_tension);
            particle_i.velo.set(
                velo_i
                    + repulsion * params.repulsion_force
                    + (average_velo - velo_i) * params.viscosity,
            );
            particle_j.velo.set(
                velo_j - repulsion * params.repulsion_force
                    + (average_velo - velo_j) * params.viscosity,
            );
        }
    }

    fn update_particles(&mut self) {
        match self.neighbor_mode {
            NeighborMode::BruteForce => self.neighbor_payload = NeighborPayload::BruteForce,
            NeighborMode::HashMap => {
                if !matches!(self.neighbor_payload, NeighborPayload::HashMap(_)) {
                    // We keep the hash map although we rebuild it every frame in the hope that some of the heap memory can be reused.
                    self.neighbor_payload = NeighborPayload::HashMap(HashMap::new());
                }
            }
            NeighborMode::SortMap => {
                if !matches!(self.neighbor_payload, NeighborPayload::SortMap { .. }) {
                    self.neighbor_payload = NeighborPayload::SortMap {
                        hash_table: vec![HashEntry::default(); self.particles.len()],
                        start_offsets: vec![usize::MAX; self.particles.len()],
                    };
                }
            }
        }

        match self.neighbor_payload {
            NeighborPayload::BruteForce => {
                let (_, time) = crate::measure_time(|| {
                    for (i, particle_i) in self.particles.iter().enumerate() {
                        for (j, particle_j) in self.particles.iter().enumerate() {
                            if i == j {
                                continue;
                            }
                            Self::update_speed(particle_i, particle_j, &self.params);
                        }
                    }
                });
                self.time_history.push_back(time);
                if MAX_HISTORY < self.time_history.len() {
                    self.time_history.pop_front();
                }
            }
            NeighborPayload::HashMap(ref mut hash_map) => {
                let (_, time) = crate::measure_time(|| {
                    hash_map.clear();
                    for (i, particle_i) in self.particles.iter().enumerate() {
                        let pos = particle_i.pos.get();
                        let grid_pos = (
                            pos.x.div_euclid(CELL_SIZE) as i32,
                            pos.y.div_euclid(CELL_SIZE) as i32,
                        );
                        hash_map.entry(grid_pos).or_default().push(i);
                    }
                });
                self.hash_time_history.push_back(time);
                if MAX_HISTORY < self.hash_time_history.len() {
                    self.hash_time_history.pop_front();
                }

                for (i, particle_i) in self.particles.iter().enumerate() {
                    let pos = particle_i.pos.get();
                    let grid_pos = (
                        pos.x.div_euclid(CELL_SIZE) as i32,
                        pos.y.div_euclid(CELL_SIZE) as i32,
                    );
                    for cy in (grid_pos.1 - 1)..=(grid_pos.1 + 1) {
                        for cx in (grid_pos.0 - 1)..=(grid_pos.0 + 1) {
                            let Some(cell) = hash_map.get(&(cx, cy)) else {
                                continue;
                            };
                            for &j in cell {
                                if i == j {
                                    continue;
                                }
                                let particle_j = &self.particles[j];
                                Self::update_speed(particle_i, particle_j, &self.params);
                            }
                        }
                    }
                }
            }
            NeighborPayload::SortMap {
                ref mut hash_table,
                ref mut start_offsets,
            } => {
                let hasher = |grid_pos: (i32, i32)| -> usize {
                    (grid_pos.0 + grid_pos.1 * 32121).rem_euclid(self.particles.len() as i32)
                        as usize
                };
                let (_, hash_time) = measure_time(|| {
                    for (i, particle_i) in self.particles.iter().enumerate() {
                        let pos = particle_i.pos.get();
                        let grid_pos = (
                            pos.x.div_euclid(CELL_SIZE) as i32,
                            pos.y.div_euclid(CELL_SIZE) as i32,
                        );
                        // let mut hasher = std::hash::SipHasher::new();
                        // grid_pos.hash(&mut hasher);
                        // let grid_hash: usize = hasher.into() % self.particle_hash.len();
                        let cell_hash = hasher(grid_pos);
                        hash_table[i] = HashEntry {
                            particle_idx: i,
                            cell_hash,
                        };
                    }
                });
                let (_, sort_time) = crate::measure_time(|| {
                    hash_table.sort_unstable_by_key(|entry| entry.cell_hash);
                    for i in 0..self.particles.len() {
                        let first = hash_table
                            .iter()
                            .enumerate()
                            .find(|(_, p)| p.cell_hash == i);
                        if let Some(first) = first {
                            start_offsets[i] = first.0;
                        }
                    }
                });
                self.sort_time_history.push_back((hash_time, sort_time));
                if MAX_HISTORY < self.sort_time_history.len() {
                    self.sort_time_history.pop_front();
                }

                for (i, particle_i) in self.particles.iter().enumerate() {
                    let pos = particle_i.pos.get();
                    let grid_pos = (
                        pos.x.div_euclid(CELL_SIZE) as i32,
                        pos.y.div_euclid(CELL_SIZE) as i32,
                    );
                    for cy in (grid_pos.1 - 1)..=(grid_pos.1 + 1) {
                        for cx in (grid_pos.0 - 1)..=(grid_pos.0 + 1) {
                            let cell_hash = hasher((cx, cy));
                            let Some(&cell_start) = start_offsets.get(cell_hash) else {
                                continue;
                            };
                            if cell_start == usize::MAX {
                                continue;
                            }
                            for entry in &hash_table[cell_start..] {
                                if i == entry.particle_idx {
                                    continue;
                                }
                                if entry.cell_hash != cell_hash {
                                    break;
                                }
                                let particle_j = &self.particles[entry.particle_idx];
                                Self::update_speed(particle_i, particle_j, &self.params);
                            }
                        }
                    }
                }

                // static INIT: AtomicBool = AtomicBool::new(false);
                // if !INIT.load(std::sync::atomic::Ordering::Relaxed) {
                //     println!("particle_hash:");
                //     for h in &self.particle_hash {
                //         println!("  {}", h);
                //     }
                //     INIT.store(true, std::sync::atomic::Ordering::Relaxed);
                //     println!("start_offsets:");
                //     for so in &self.start_offsets {
                //         println!("  {}", so);
                //     }
                //     INIT.store(true, std::sync::atomic::Ordering::Relaxed);
                // }
            }
        }

        if let Some((mouse_pos, delta)) = self.mouse_down {
            for particle_i in self.particles.iter() {
                let pos_i = particle_i.pos.get();
                let dist2 = (pos_i - mouse_pos.to_vec2()).length_sq();
                if dist2 < self.mouse_effect_radius.powi(2) {
                    particle_i.velo.set(delta);
                }
            }
        }

        for particle in &self.particles {
            let pos = particle.pos.get();
            let mut velo = particle.velo.get();
            velo.y -= self.gravity;
            let newpos = pos + velo;
            let croppos = pos2(
                newpos.x.min(self.rect.max.x).max(self.rect.min.x),
                newpos.y.min(self.rect.max.y).max(self.rect.min.y),
            );
            particle.pos.set(croppos.to_vec2());
            if newpos.x < self.rect.min.x && velo.x < 0. {
                velo.x = -velo.x * self.restitution;
            }
            if self.rect.max.x < newpos.x && 0. < velo.x {
                velo.x = -velo.x * self.restitution;
            }
            if newpos.y < self.rect.min.y && velo.y < 0. {
                velo.y = -velo.y * self.restitution;
            }
            if self.rect.max.y < newpos.y && 0. < velo.y {
                velo.y = -velo.y * self.restitution;
            }
            particle.velo.set(velo);
        }
    }

    fn ui_panel(&mut self, ui: &mut Ui) {
        if ui.button("Reset").clicked() {
            self.reset();
        }
        ui.label("Num particles (needs reset):");
        ui.add(egui::widgets::Slider::new(
            &mut self.num_particles,
            1..=15000,
        ));
        ui.label("Width:");
        if ui
            .add(egui::widgets::Slider::new(
                &mut self.rect.max.x,
                1.0..=100.0,
            ))
            .changed()
        {
            self.rect.min.x = -self.rect.max.x;
        };
        ui.label("Height:");
        if ui
            .add(egui::widgets::Slider::new(
                &mut self.rect.max.y,
                1.0..=100.0,
            ))
            .changed()
        {
            self.rect.min.y = -self.rect.max.y;
        };
        ui.label("Restitution:");
        ui.add(egui::widgets::Slider::new(&mut self.restitution, (0.)..=1.));
        ui.label("Repulsion force:");
        ui.add(egui::widgets::Slider::new(
            &mut self.params.repulsion_force,
            (0.)..=1.0,
        ));
        ui.label("Viscosity:");
        ui.add(egui::widgets::Slider::new(
            &mut self.params.viscosity,
            (0.)..=0.01,
        ));
        ui.label("Surface tension:");
        ui.add(egui::widgets::Slider::new(
            &mut self.params.surface_tension,
            (0.)..=0.5,
        ));
        ui.label("Gravity:");
        ui.add(egui::widgets::Slider::new(&mut self.gravity, (0.)..=0.1));
        ui.label("Mouse effect radius:");
        ui.add(egui::widgets::Slider::new(
            &mut self.mouse_effect_radius,
            (0.)..=10.,
        ));
        ui.group(|ui| {
            ui.label("Neighbor search:");
            ui.radio_value(
                &mut self.neighbor_mode,
                NeighborMode::BruteForce,
                "Brute force",
            );
            ui.radio_value(&mut self.neighbor_mode, NeighborMode::HashMap, "Hash map");
            ui.radio_value(&mut self.neighbor_mode, NeighborMode::SortMap, "Sort map");
        });
        ui.checkbox(&mut self.show_grid, "Show grid");
        ui.checkbox(&mut self.show_grid_count, "Show grid count");
        ui.checkbox(&mut self.color_by_speed, "Color by speed");
        ui.checkbox(&mut self.show_time_plot, "Show time plot");
    }

    fn plot_time(&mut self, ui: &mut Ui) {
        let plot = Plot::new("plot");
        plot.legend(Legend::default()).show(ui, |plot_ui| {
            let gen_line = |points, i, name| {
                Line::new(points)
                    .color(eframe::egui::Color32::from_rgb(
                        (i % 2 * 200) as u8,
                        (i % 4 * 200) as u8,
                        (i % 8 * 100) as u8,
                    ))
                    .name(name)
            };
            let points: PlotPoints = self
                .time_history
                .iter()
                .enumerate()
                .map(|(t, v)| [t as f64, *v])
                .collect();
            plot_ui.line(gen_line(points, 0, "BruteForce time"));
            let points: PlotPoints = self
                .hash_time_history
                .iter()
                .enumerate()
                .map(|(t, v)| [t as f64, *v])
                .collect();
            plot_ui.line(gen_line(points, 1, "HashMap time"));
            let hash_points: PlotPoints = self
                .sort_time_history
                .iter()
                .enumerate()
                .map(|(t, v)| [t as f64, v.0])
                .collect();
            let sort_points: PlotPoints = self
                .sort_time_history
                .iter()
                .enumerate()
                .map(|(t, v)| [t as f64, v.1])
                .collect();
            plot_ui.line(gen_line(hash_points, 2, "SortMap Hash time"));
            plot_ui.line(gen_line(sort_points, 3, "SortMap sort time"));
        });
    }
}

impl eframe::App for RusHydroApp {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        ctx.request_repaint();

        self.update_particles();

        egui::SidePanel::right("side_panel")
            .min_width(200.)
            .show(ctx, |ui| self.ui_panel(ui));

        if self.show_time_plot {
            egui::TopBottomPanel::bottom("bottom_chart")
                .resizable(true)
                .show(ctx, |ui| {
                    self.plot_time(ui);
                });
        }

        egui::CentralPanel::default()
            // .resizable(true)
            // .min_height(100.)
            .show(ctx, |ui| {
                self.paint_canvas(ui);
            });
    }
}
