use std::cell::Cell;

use eframe::{
    egui::{self, Context, Frame, Ui},
    epaint::{pos2, vec2, Color32, Pos2, Rect, Vec2},
};
use rand::{thread_rng, Rng};

const SCALE: f32 = 10.;
const NUM_PARTICLES: usize = 200;
const PARTICLE_RADIUS: f32 = 2.;
const PARTICLE_RADIUS2: f32 = PARTICLE_RADIUS * PARTICLE_RADIUS;
const RESTITUTION: f32 = 0.5;
const REPULSION_FORCE: f32 = 0.005;
const G: f32 = 0.01;

struct Particle {
    pos: Cell<Vec2>,
    velo: Cell<Vec2>,
}

pub struct RusHydroApp {
    particles: Vec<Particle>,
    rect: Rect,
    restitution: f32,
    repulsion_force: f32,
    gravity: f32,
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
            restitution: RESTITUTION,
            repulsion_force: REPULSION_FORCE,
            gravity: G,
        }
    }

    fn paint_canvas(&mut self, ui: &mut Ui) {
        Frame::canvas(ui.style()).show(ui, |ui| {
            let (response, painter) =
                ui.allocate_painter(ui.available_size(), egui::Sense::click());

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

            for particle in &self.particles {
                painter.circle(
                    to_pos2(particle.pos.get().to_pos2()),
                    5.,
                    Color32::BLUE,
                    (1., Color32::BLACK),
                );
            }
        });
    }

    fn reset(&mut self) {
        let mut rng = thread_rng();
        let particles = (0..NUM_PARTICLES)
            .map(|_| Particle {
                pos: Cell::new(vec2(
                    rng.gen_range(self.rect.min.x..self.rect.max.x),
                    rng.gen_range(self.rect.min.y..self.rect.max.y),
                )),
                velo: Cell::new(Vec2::ZERO),
            })
            .collect();
        self.particles = particles;
    }

    fn update_particles(&mut self) {
        for (i, particle_i) in self.particles.iter().enumerate() {
            for (j, particle_j) in self.particles.iter().enumerate() {
                if i == j {
                    continue;
                }
                let pos_i = particle_i.pos.get();
                let pos_j = particle_j.pos.get();
                let delta = pos_i - pos_j;
                let dist2 = delta.length_sq();
                if 0. < dist2 && dist2 < PARTICLE_RADIUS2 {
                    let repulsing = delta / dist2.sqrt();
                    let velo_i = particle_i.velo.get();
                    particle_i
                        .velo
                        .set(velo_i + repulsing * self.repulsion_force);
                    let velo_j = particle_j.velo.get();
                    particle_j
                        .velo
                        .set(velo_j - repulsing * self.repulsion_force);
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
}

impl eframe::App for RusHydroApp {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        ctx.request_repaint();

        self.update_particles();

        eframe::egui::SidePanel::right("side_panel")
            .min_width(200.)
            .show(ctx, |ui| {
                if ui.button("Reset").clicked() {
                    self.reset();
                }
                ui.label("Restitution:");
                ui.add(egui::widgets::Slider::new(&mut self.restitution, (0.)..=1.));
                ui.label("Repulsion force:");
                ui.add(egui::widgets::Slider::new(
                    &mut self.repulsion_force,
                    (0.)..=0.01,
                ));
                ui.label("Gravity:");
                ui.add(egui::widgets::Slider::new(&mut self.gravity, (0.)..=0.1));
            });

        egui::CentralPanel::default()
            // .resizable(true)
            // .min_height(100.)
            .show(ctx, |ui| {
                self.paint_canvas(ui);
            });
    }
}
