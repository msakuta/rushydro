use eframe::epaint::{pos2, vec2, Pos2, Rect, Vec2};

use super::{RusHydroApp, PARTICLE_RENDER_RADIUS, SCALE};

#[derive(PartialEq, Eq, Clone, Copy)]
pub(crate) enum Obstacles {
    None,
    Snake,
    Slope,
}

pub(super) struct Obstacle {
    transform: Transform,
    halfsize: Vec2,
}

impl Obstacle {
    pub fn new(offset: Vec2, rotation: f32, halfsize: Vec2) -> Self {
        Self {
            transform: Transform::new(offset, rotation),
            halfsize,
        }
    }

    pub fn _offset(&self) -> Vec2 {
        self.transform.offset
    }

    pub fn _halfsize(&self) -> Vec2 {
        self.halfsize
    }

    pub fn vertices(&self) -> Vec<Pos2> {
        let min = -self.halfsize;
        let max = self.halfsize;
        [min, vec2(max.x, min.y), max, vec2(min.x, max.y)]
            .into_iter()
            .map(|v| self.transform.apply(v.to_pos2()))
            .collect()
    }

    pub fn shrink(&self, ofs: f32) -> Self {
        Self {
            transform: self.transform,
            halfsize: self.halfsize - Vec2::splat(ofs),
        }
    }

    /// Returns the distance (negative) and the normal vector of that surface.
    pub fn hit_test(&self, pos: Vec2) -> Option<(f32, Vec2)> {
        let edge_distance = |p: Vec2, org: Vec2, normal: Vec2| {
            let dp = p - org;
            dp.dot(normal)
        };
        let up = self.transform.apply_vec(vec2(0., 1.));
        let top = edge_distance(pos, self.transform.offset + up * self.halfsize.y, up);
        let down = -up;
        let bottom = edge_distance(pos, self.transform.offset + down * self.halfsize.y, down);
        let left_dir = self.transform.apply_vec(vec2(1., 0.));
        let left = edge_distance(
            pos,
            self.transform.offset + left_dir * self.halfsize.x,
            left_dir,
        );
        let right_dir = -left_dir;
        let right = edge_distance(
            pos,
            self.transform.offset + right_dir * self.halfsize.x,
            right_dir,
        );

        if 0. < left || 0. < top || 0. < right || 0. < bottom {
            return None;
        }

        [
            (left, left_dir),
            (top, up),
            (right, right_dir),
            (bottom, down),
        ]
        .iter()
        .enumerate()
        .fold(None, |acc: Option<(usize, (f32, Vec2))>, cur| {
            if 0. < cur.1 .0 {
                return acc;
            }
            if let Some(acc) = acc {
                if acc.1 .0 < cur.1 .0 {
                    Some((cur.0, *cur.1))
                } else {
                    Some(acc)
                }
            } else {
                Some((cur.0, *cur.1))
            }
        })
        .map(|(_, v)| v)
    }
}

/// Rotate and offset. No scaling.
#[derive(Clone, Copy, Debug)]
struct Transform {
    offset: Vec2,
    rotation: f32,
}

impl Transform {
    fn new(offset: Vec2, rotation: f32) -> Self {
        Self { offset, rotation }
    }

    fn apply(&self, pos: Pos2) -> Pos2 {
        let s = self.rotation.sin();
        let c = self.rotation.cos();
        pos2(
            pos.x * c + pos.y * s + self.offset.x,
            pos.x * -s + pos.y * c + self.offset.y,
        )
    }

    fn apply_vec(&self, v: Vec2) -> Vec2 {
        let s = self.rotation.sin();
        let c = self.rotation.cos();
        vec2(v.x * c + v.y * s, v.x * -s + v.y * c)
    }
}

impl RusHydroApp {
    pub(super) const SLOPE_ANGLE: f32 = std::f32::consts::PI * 0.075;

    pub(super) fn gen_obstacles(rect: &Rect, obs: Obstacles) -> Vec<Obstacle> {
        match obs {
            Obstacles::Snake => {
                const OFFSET: f32 = PARTICLE_RENDER_RADIUS / SCALE;
                vec![
                    Obstacle::new(
                        vec2(-rect.width() * 0.25 - OFFSET, rect.height() / 4.),
                        0.,
                        vec2(rect.width() * 0.25, 2.),
                    ),
                    Obstacle::new(
                        vec2(rect.width() * 0.25 + OFFSET, -rect.height() / 4.),
                        0.,
                        vec2(rect.width() * 0.25, 2.),
                    ),
                ]
            }
            Obstacles::Slope => {
                let y_off = rect.width() * 0.5 * Self::SLOPE_ANGLE.sin();
                vec![Obstacle::new(
                    vec2(0., -rect.height() + y_off),
                    Self::SLOPE_ANGLE,
                    vec2(rect.width(), rect.height() * 0.5),
                )]
            }
            _ => vec![],
        }
    }
}
