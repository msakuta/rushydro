use eframe::epaint::{pos2, vec2, Pos2, Rect, Vec2};

use super::{RusHydroApp, PARTICLE_RENDER_RADIUS, SCALE};

#[derive(PartialEq, Eq, Clone, Copy)]
pub(crate) enum Obstacles {
    None,
    Snake,
    Slope,
    WaterMill,
}

pub(super) struct RectObstacle {
    transform: Transform,
    halfsize: Vec2,
    velo: Vec2,
    angular_velo: f32,
    mass: f32,
    /// Moment of Inertia, kg * m^2
    moi: f32,
    pivot: Option<Vec2>,
    sense_gravity: bool,
}

impl RectObstacle {
    pub fn new(offset: Vec2, rotation: f32, halfsize: Vec2, mass: f32) -> Self {
        Self {
            transform: Transform::new(offset, rotation),
            halfsize,
            velo: Vec2::ZERO,
            angular_velo: 0.,
            mass,
            moi: mass / 12. * (halfsize.x.powi(2) + halfsize.y.powi(2)),
            pivot: None,
            sense_gravity: false,
        }
    }

    pub fn with_pivot(self, pivot: Vec2) -> Self {
        Self {
            pivot: Some(pivot),
            ..self
        }
    }

    pub fn _with_gravity(self) -> Self {
        Self {
            sense_gravity: true,
            ..self
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
            halfsize: self.halfsize - Vec2::splat(ofs),
            ..*self
        }
    }

    pub fn update(&mut self, gravity: f32, delta_time: f32) {
        if self.mass != 0. {
            if self.sense_gravity {
                self.velo.y -= gravity * delta_time;
            }
            self.transform.offset += self.velo * delta_time;
            self.transform.rotation += self.angular_velo * delta_time;
        }
    }
}

pub(super) struct ObstacleHitResult {
    pub dist: f32,
    pub normal: Vec2,
    pub velo: Vec2,
}

impl RectObstacle {
    /// Returns the distance (negative) and the normal vector of that surface.
    pub fn hit_test(&self, pos: Vec2) -> Option<ObstacleHitResult> {
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
        .map(|(_, v)| {
            let dp = pos - self.transform.offset;
            ObstacleHitResult {
                dist: v.0,
                normal: v.1,
                velo: self.angular_velo * rotate90(dp),
            }
        })
    }

    pub(super) fn apply_impulse(&mut self, impulse: Vec2, pos: Vec2) {
        if self.mass != 0. {
            if self.pivot.is_none() {
                self.velo += impulse / self.mass;
            } else {
                let moment = cross(pos - self.transform.offset, impulse);
                self.angular_velo += moment / self.moi;
            }
        }
    }
}

fn rotate90(v: Vec2) -> Vec2 {
    Vec2 { x: v.y, y: -v.x }
}

fn cross(a: Vec2, b: Vec2) -> f32 {
    a.y * b.x - a.x * b.y
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
    pub(super) const SLOPE_ANGLE: f32 = std::f32::consts::PI * 0.05;

    pub(super) fn gen_obstacles(rect: &Rect, obs: Obstacles) -> Vec<RectObstacle> {
        let gen_slope = || {
            let y_off = rect.width() * 0.5 * Self::SLOPE_ANGLE.sin();
            RectObstacle::new(
                vec2(0., -rect.height() + y_off),
                Self::SLOPE_ANGLE,
                vec2(rect.width(), rect.height() * 0.5),
                0.,
            )
        };
        match obs {
            Obstacles::Snake => {
                const OFFSET: f32 = PARTICLE_RENDER_RADIUS / SCALE;
                vec![
                    RectObstacle::new(
                        vec2(-rect.width() * 0.25 - OFFSET, rect.height() / 4.),
                        Self::SLOPE_ANGLE,
                        vec2(rect.width() * 0.3, 2.),
                        0.,
                    ),
                    RectObstacle::new(
                        vec2(rect.width() * 0.25 + OFFSET, -rect.height() / 4.),
                        -Self::SLOPE_ANGLE,
                        vec2(rect.width() * 0.3, 2.),
                        0.,
                    ),
                ]
            }
            Obstacles::Slope => {
                vec![gen_slope()]
            }
            Obstacles::WaterMill => {
                let axis = vec2(0., rect.height() * 0.1);
                vec![
                    gen_slope(),
                    RectObstacle::new(axis, 0., vec2(2., rect.height() * 0.4), 1e4)
                        .with_pivot(Vec2::ZERO),
                ]
            }
            _ => vec![],
        }
    }

    pub(super) fn update_obstacles(&mut self, delta_time: f32) {
        for obstacle in &mut self.obstacles {
            obstacle.update(self.gravity, delta_time);
        }
    }
}
