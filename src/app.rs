mod particles;

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

use self::particles::Particle;

use crate::marching_squares::{border_pixel, cell_border_index, pick_bits, Shape};

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
    show_particles: bool,
    show_surface: bool,
    show_grid: bool,
    show_grid_count: bool,
    color_by_speed: bool,
    show_time_plot: bool,
    time_history: VecDeque<f64>,
    hash_time_history: VecDeque<f64>,
    sort_time_history: VecDeque<(f64, f64)>,
    density_map: Vec<f32>,
    density_shape: Shape,
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
        let (density_map, density_shape) = Self::gen_board(&rect);
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
            show_particles: true,
            show_surface: true,
            show_grid: false,
            show_grid_count: false,
            color_by_speed: true,
            show_time_plot: false,
            time_history: VecDeque::new(),
            hash_time_history: VecDeque::new(),
            sort_time_history: VecDeque::new(),
            density_map,
            density_shape,
        }
    }

    fn gen_board(rect: &Rect) -> (Vec<f32>, Shape) {
        let width = (rect.width() / CELL_SIZE) as usize + 1;
        let height = (rect.height() / CELL_SIZE) as usize + 1;
        (vec![0.; width * height], (width as isize, height as isize))
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

            if self.show_surface {
                self.density_map.fill(0.);
                for particle in self.particles.iter() {
                    let density_pos = (
                        (particle.pos.get().x - self.rect.min.x).div_euclid(CELL_SIZE) as isize,
                        (particle.pos.get().y - self.rect.min.y).div_euclid(CELL_SIZE) as isize,
                    );
                    if 0 <= density_pos.0
                        && density_pos.0 < self.density_shape.0
                        && 0 <= density_pos.1
                        && density_pos.1 < self.density_shape.1
                    {
                        self.density_map
                            [(density_pos.0 + density_pos.1 * self.density_shape.0) as usize] += 1.;
                    }
                }

                for cy in 0..self.density_shape.1 - 1 {
                    let offset_y = cy as f32 * CELL_SIZE + self.rect.min.y;
                    for cx in 0..self.density_shape.0 - 1 {
                        let offset_x = cx as f32 * CELL_SIZE + self.rect.min.x;
                        let bits = pick_bits(&self.density_map, self.density_shape, (cx, cy), 0.5);
                        if !border_pixel(bits) {
                            continue;
                        }
                        if let Some(lines) = cell_border_index(bits) {
                            for line in lines.chunks(4) {
                                let points = [
                                    to_pos2(pos2(line[0] + offset_x, line[1] + offset_y)),
                                    to_pos2(pos2(line[2] + offset_x, line[3] + offset_y)),
                                ];
                                // let poly = PathShape::convex_polygon(points.clone(), Color32::from_rgb(127, 191, 255), (0., Color32::BLUE));
                                // painter.add(poly);
                                // let stroke = PathShape::closed_line(points, );
                                painter.line_segment(points, (1., Color32::RED));
                            }
                        }
                    }
                }
            }

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

            if self.show_particles {
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
            }
        });
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
            (self.density_map, self.density_shape) = Self::gen_board(&self.rect);
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
            (self.density_map, self.density_shape) = Self::gen_board(&self.rect);
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
        ui.checkbox(&mut self.show_particles, "Show particles");
        ui.checkbox(&mut self.show_surface, "Show surface");
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
