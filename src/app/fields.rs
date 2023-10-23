use eframe::{
    egui::{self, Painter, TextureOptions},
    epaint::{pos2, Color32, PathShape, Pos2, Rect, Stroke},
};

use crate::marching_squares::{
    border_pixel, cell_border_interpolated, cell_polygon_interpolated, pick_bits, pick_values,
    Shape,
};

use super::particles::Particle;

pub(super) struct Fields {
    shape: Shape,
    resolution: f32,
    density_map: Vec<f32>,
    temperature_map: Vec<f32>,
}

impl Fields {
    pub fn new(rect: &Rect, resolution: f32) -> Self {
        let width = (rect.width() / resolution) as usize + 1;
        let height = (rect.height() / resolution) as usize + 1;
        let map = vec![0.; width * height];
        Self {
            shape: (width as isize, height as isize),
            resolution,
            density_map: map.clone(),
            temperature_map: map,
        }
    }

    pub fn update(&mut self, rect: &Rect, particles: &[Particle], radius: f32) {
        self.density_map.fill(0.);
        self.temperature_map.fill(0.);
        let resol = self.resolution;
        let pix_rad = (radius / resol).ceil() as isize;
        for particle in particles.iter() {
            let pos = particle.pos.get();
            let density_pos = ((pos.x - rect.min.x) / resol, (pos.y - rect.min.y) / resol);
            let density_idx = (
                density_pos.0.floor() as isize,
                density_pos.1.floor() as isize,
            );
            for cy in density_idx.1 - pix_rad..=density_idx.1 + pix_rad {
                for cx in density_idx.0 - pix_rad..=density_idx.0 + pix_rad {
                    if 0 <= cx && cx < self.shape.0 && 0 <= cy && cy < self.shape.1 {
                        let dx = cx as f32 - density_pos.0;
                        let dy = cy as f32 - density_pos.1;
                        let my_density = 1. / (1. + dx * dx + dy * dy);
                        let buf_idx = (cx + cy * self.shape.0) as usize;
                        let cell_temp = &mut self.temperature_map[buf_idx];
                        let cell_dens = &mut self.density_map[buf_idx];
                        *cell_temp = (*cell_temp * *cell_dens + particle.temp.get() * my_density)
                            / (*cell_dens + my_density);
                        *cell_dens += my_density;
                    }
                }
            }
        }
    }

    pub fn render_image(&self, painter: &Painter, rect: &Rect, to_pos2: &impl Fn(Pos2) -> Pos2) {
        let bits_x = self.shape.0 - 1;
        let bits_y = self.shape.1 - 1;
        let mut bits = vec![0u8; (bits_x * bits_y) as usize];
        for cy in 0..bits_y {
            for cx in 0..bits_x {
                bits[(cx + cy * bits_x) as usize] =
                    pick_bits(&self.density_map, self.shape, (cx, cy), 0.5);
            }
        }

        let image = bits
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let x = i % bits_x as usize;
                let y = i / bits_x as usize;
                if 15 == *v {
                    Self::temperature_rgb(self.temperature_map[x + y * self.shape.0 as usize])
                } else if (x + y) % 2 == 0 {
                    [255, 251, 251]
                } else {
                    [255; 3]
                }
            })
            .flatten()
            .collect::<Vec<_>>();

        let image = egui::ColorImage::from_rgb([bits_x as usize, bits_y as usize], &image);

        let texture = painter.ctx().load_texture(
            "my-image",
            image,
            TextureOptions {
                magnification: egui::TextureFilter::Nearest,
                minification: egui::TextureFilter::Linear,
            },
        );

        const UV: Rect = Rect::from_min_max(pos2(0., 1.), Pos2::new(1.0, 0.0));
        let mut scr_rect = Rect::from_min_max(
            to_pos2(rect.min),
            to_pos2(pos2(
                rect.min.x + bits_x as f32 * self.resolution,
                rect.min.y + bits_y as f32 * self.resolution,
            )),
        );
        std::mem::swap(&mut scr_rect.min.y, &mut scr_rect.max.y);
        painter.image(texture.id(), scr_rect, UV, Color32::WHITE);
    }

    pub fn render_surface(
        &self,
        painter: &Painter,
        rect: &Rect,
        fill_polygon: bool,
        to_pos2: &impl Fn(Pos2) -> Pos2,
    ) {
        let resol = self.resolution;
        for cy in 0..self.shape.1 - 1 {
            let offset_y = (cy as f32 + 0.5) * resol + rect.min.y;
            for cx in 0..self.shape.0 - 1 {
                let offset_x = (cx as f32 + 0.5) * resol + rect.min.x;
                let bits = pick_bits(&self.density_map, self.shape, (cx, cy), 0.5);
                if !border_pixel(bits) {
                    continue;
                }
                let values = pick_values(&self.density_map, self.shape, (cx, cy));
                if fill_polygon {
                    cell_polygon_interpolated(bits, values, |points| {
                        let points = points
                            .chunks(2)
                            .map(|x| {
                                to_pos2(pos2(
                                    x[0] * resol * 0.5 + offset_x,
                                    x[1] * resol * 0.5 + offset_y,
                                ))
                            })
                            .collect();
                        let rgb = Self::temperature_rgb(
                            self.temperature_map[(cx + cy * self.shape.0) as usize],
                        );
                        let poly = PathShape::convex_polygon(
                            points,
                            Color32::from_rgb(rgb[0], rgb[1], rgb[2]),
                            Stroke::NONE,
                        );
                        painter.add(poly);
                    });
                }
                if let Some((lines, len)) = cell_border_interpolated(bits, values) {
                    for line in lines.chunks(4).take(len / 4) {
                        let points = [
                            to_pos2(pos2(
                                line[0] * resol * 0.5 + offset_x,
                                line[1] * resol * 0.5 + offset_y,
                            )),
                            to_pos2(pos2(
                                line[2] * resol * 0.5 + offset_x,
                                line[3] * resol * 0.5 + offset_y,
                            )),
                        ];
                        painter.line_segment(points, (1., Color32::BLUE));
                    }
                }
            }
        }
    }

    fn temperature_rgb(temperature: f32) -> [u8; 3] {
        let red = (255. * (temperature * 0.5 + 0.5)).clamp(0., 255.) as u8;
        let gb = (255. * (1. - temperature * 0.5)).clamp(0., 255.) as u8;
        [red, gb, gb]
    }
}
