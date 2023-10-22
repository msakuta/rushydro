use eframe::epaint::{pos2, vec2, Pos2, Rect, Vec2};

use super::{RusHydroApp, PARTICLE_RENDER_RADIUS, SCALE};

#[derive(PartialEq, Eq, Clone, Copy)]
pub(crate) enum Obstacles {
    None,
    Snake,
    Slope,
    WaterMill,
}

pub(super) struct Obstacle {
    shape: ObstacleShape,
    transform: Transform,
    velo: Vec2,
    angular_velo: f32,
    mass: f32,
    /// Moment of Inertia, kg * m^2
    moi: f32,
    pivot: Option<Vec2>,
    sense_gravity: bool,
    hydrophilic: bool,
}

pub(super) enum ObstacleShape {
    Rect(RectObstacle),
    Compound(Vec<ObstacleShape>),
}

impl Obstacle {
    pub fn new(ro: RectObstacle) -> Self {
        let mass = ro.mass;
        let moi = ro.moi();
        Self {
            shape: ObstacleShape::Rect(ro),
            transform: Transform::default(),
            velo: Vec2::ZERO,
            angular_velo: 0.,
            mass,
            moi,
            pivot: None,
            sense_gravity: false,
            hydrophilic: false,
        }
    }

    pub fn new_compound(v: Vec<ObstacleShape>, pos: Vec2) -> Self {
        let mass = v.iter().map(ObstacleShape::mass).sum();
        let moi = v
            .iter()
            .map(|c| {
                // TODO: fix moment of inertia around non-centroid
                let local_moi = c.moi();
                local_moi
            })
            .sum();
        Self {
            shape: ObstacleShape::Compound(v),
            transform: Transform {
                offset: pos,
                rotation: 0.,
            },
            velo: Vec2::ZERO,
            angular_velo: 0.,
            mass,
            moi,
            pivot: None,
            sense_gravity: false,
            hydrophilic: false,
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

    pub fn with_hydrophilic(self, v: bool) -> Self {
        Self {
            hydrophilic: v,
            ..self
        }
    }

    pub fn is_hydrophilic(&self) -> bool {
        self.hydrophilic
    }

    /// Call the callback for each of the shape in a potentially compound object.
    ///
    /// You can pass offset to render the shape with particle size into account.
    pub fn foreach_shape(&self, f: &mut impl FnMut(&[Pos2]), ofs: f32) {
        self.shape.foreach_shape(
            &mut |points| {
                let points: Vec<_> = points
                    .iter()
                    .map(|pos| self.transform.apply(*pos))
                    .collect();
                f(&points)
            },
            ofs,
        );
    }

    /// Returns the distance (negative) and the normal vector of that surface.
    ///
    /// hit_test is similar to `distance`, but the difference is that hit test only returns if the given point
    /// intersects with the shape, thus utilizing some optimizaion to make the calculation faster.
    pub fn hit_test(&self, pos: Vec2) -> Option<ObstacleHitResult> {
        let local_pos = self.transform.apply_inverse(pos.to_pos2());
        self.shape
            .distance(local_pos.to_vec2(), true)
            .map(|mut res| {
                res.normal = self.transform.apply_vec(res.normal);
                let velo = self.angular_velo * rotate90(local_pos.to_vec2());
                res.velo = self.transform.apply_vec(velo);
                res
            })
    }

    /// Returns the distance (either positive or negative) and the normal vector of that surface.
    pub fn distance(&self, pos: Vec2) -> Option<ObstacleHitResult> {
        let local_pos = self.transform.apply_inverse(pos.to_pos2());
        self.shape
            .distance(local_pos.to_vec2(), false)
            .map(|mut res| {
                res.normal = self.transform.apply_vec(res.normal);
                let velo = self.angular_velo * rotate90(local_pos.to_vec2());
                res.velo = self.transform.apply_vec(velo);
                res
            })
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

impl ObstacleShape {
    pub fn mass(&self) -> f32 {
        match self {
            Self::Rect(r) => r.mass,
            Self::Compound(c) => c.iter().map(|c| c.mass()).sum(),
        }
    }

    fn moi(&self) -> f32 {
        match self {
            Self::Rect(r) => r.moi(),
            Self::Compound(c) => c.iter().map(|c| c.moi()).sum(),
        }
    }

    fn foreach_shape(&self, f: &mut impl FnMut(&[Pos2]), ofs: f32) {
        match self {
            Self::Rect(r) => f(&r.vertices(ofs)),
            Self::Compound(c) => {
                for v in c {
                    v.foreach_shape(f, ofs);
                }
            }
        }
    }

    fn distance(&self, pos: Vec2, only_hits: bool) -> Option<ObstacleHitResult> {
        match self {
            Self::Rect(r) => r.hit_test(pos, only_hits),
            Self::Compound(c) => c.iter().fold(None, |acc, cur| {
                let Some(res) = cur.distance(pos, only_hits) else {
                    return acc;
                };
                if let Some(acc) = acc {
                    if res.dist < acc.dist {
                        Some(res)
                    } else {
                        Some(acc)
                    }
                } else {
                    Some(res)
                }
            }),
        }
    }
}

pub(super) struct RectObstacle {
    transform: Transform,
    halfsize: Vec2,
    mass: f32,
}

impl RectObstacle {
    pub fn new(offset: Vec2, rotation: f32, halfsize: Vec2, mass: f32) -> Self {
        Self {
            transform: Transform::new(offset, rotation),
            halfsize,
            mass,
        }
    }

    pub fn _offset(&self) -> Vec2 {
        self.transform.offset
    }

    pub fn _halfsize(&self) -> Vec2 {
        self.halfsize
    }

    /// Moment of Inertia around center of gravity.
    fn moi(&self) -> f32 {
        self.mass / 12. * (self.halfsize.x.powi(2) + self.halfsize.y.powi(2))
    }

    pub fn vertices(&self, ofs: f32) -> Vec<Pos2> {
        let min = -self.halfsize - Vec2::splat(ofs);
        let max = self.halfsize + Vec2::splat(ofs);
        [min, vec2(max.x, min.y), max, vec2(min.x, max.y)]
            .into_iter()
            .map(|v| self.transform.apply(v.to_pos2()))
            .collect()
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct ObstacleHitResult {
    pub dist: f32,
    pub normal: Vec2,
    pub velo: Vec2,
}

impl PartialEq for ObstacleHitResult {
    fn eq(&self, other: &Self) -> bool {
        self.dist == other.dist
    }
}

impl PartialOrd for ObstacleHitResult {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.dist.partial_cmp(&other.dist)
    }
}

impl RectObstacle {
    /// Returns the distance (negative) and the normal vector of that surface.
    pub fn hit_test(&self, pos: Vec2, only_hits: bool) -> Option<ObstacleHitResult> {
        let get_edge_distance = |edge_dist: f32, local_normal: Vec2| {
            let normal = self.transform.apply_vec(local_normal);
            let org = self.transform.offset + normal * edge_dist;
            let dp = pos - org;
            (dp.dot(normal), normal)
        };
        let top = get_edge_distance(self.halfsize.y, vec2(0., 1.));
        let bottom = get_edge_distance(self.halfsize.y, vec2(0., -1.));
        let left = get_edge_distance(self.halfsize.x, vec2(1., 0.));
        let right = get_edge_distance(self.halfsize.x, vec2(-1., 0.));

        if only_hits && (0. < left.0 || 0. < top.0 || 0. < right.0 || 0. < bottom.0) {
            return None;
        }

        if !only_hits {
            let local_pos = self.transform.apply_inverse(pos.to_pos2());
            if self.halfsize.x < local_pos.x.abs() && self.halfsize.y < local_pos.y.abs() {
                let corner_pos = pos2(
                    local_pos.x.signum() * self.halfsize.x,
                    local_pos.y.signum() * self.halfsize.y,
                );
                let delta = self.transform.apply_vec(local_pos - corner_pos);
                let len = delta.length();
                return Some(ObstacleHitResult {
                    dist: len,
                    normal: delta / len,
                    velo: Vec2::ZERO,
                });
            }
        }

        [left, top, right, bottom]
            .iter()
            .fold(None, |acc: Option<(f32, Vec2)>, cur: &(f32, Vec2)| {
                if only_hits && 0. < cur.0 {
                    return acc;
                }
                if let Some(acc) = acc {
                    if acc.0 < cur.0 {
                        Some(*cur)
                    } else {
                        Some(acc)
                    }
                } else {
                    Some(*cur)
                }
            })
            .map(|v| ObstacleHitResult {
                dist: v.0,
                normal: v.1,
                velo: Vec2::ZERO,
            })
    }
}

fn rotate90(v: Vec2) -> Vec2 {
    Vec2 { x: v.y, y: -v.x }
}

fn cross(a: Vec2, b: Vec2) -> f32 {
    a.y * b.x - a.x * b.y
}

/// Rotate and offset. No scaling.
#[derive(Clone, Copy, Debug, Default)]
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

    fn apply_inverse(&self, mut pos: Pos2) -> Pos2 {
        let s = self.rotation.sin();
        let c = self.rotation.cos();
        pos.x -= self.offset.x;
        pos.y -= self.offset.y;
        pos2(pos.x * c + pos.y * -s, pos.x * s + pos.y * c)
    }
}

impl RusHydroApp {
    pub(super) const SLOPE_ANGLE: f32 = std::f32::consts::PI * 0.05;

    pub(super) fn gen_obstacles(rect: &Rect, obs: Obstacles) -> Vec<Obstacle> {
        let gen_slope = || {
            let y_off = rect.width() * 0.5 * Self::SLOPE_ANGLE.sin();
            Obstacle::new(RectObstacle::new(
                vec2(0., -rect.height() + y_off),
                Self::SLOPE_ANGLE,
                vec2(rect.width(), rect.height() * 0.5),
                0.,
            ))
        };
        match obs {
            Obstacles::Snake => {
                const OFFSET: f32 = PARTICLE_RENDER_RADIUS / SCALE;
                vec![
                    Obstacle::new(RectObstacle::new(
                        vec2(-rect.width() * 0.25 - OFFSET, rect.height() / 4.),
                        Self::SLOPE_ANGLE,
                        vec2(rect.width() * 0.3, 2.),
                        0.,
                    )),
                    Obstacle::new(RectObstacle::new(
                        vec2(rect.width() * 0.25 + OFFSET, -rect.height() / 4.),
                        -Self::SLOPE_ANGLE,
                        vec2(rect.width() * 0.3, 2.),
                        0.,
                    )),
                ]
            }
            Obstacles::Slope => {
                vec![gen_slope()]
            }
            Obstacles::WaterMill => {
                let axis = vec2(0., rect.height() * 0.1);
                vec![
                    gen_slope(),
                    Obstacle::new_compound(
                        vec![
                            ObstacleShape::Rect(RectObstacle::new(
                                Vec2::ZERO,
                                0.,
                                vec2(2., rect.height() * 0.4),
                                1e4,
                            )),
                            ObstacleShape::Rect(RectObstacle::new(
                                Vec2::ZERO,
                                std::f32::consts::PI * 0.5,
                                vec2(2., rect.height() * 0.4),
                                1e4,
                            )),
                        ],
                        axis,
                    )
                    .with_pivot(Vec2::ZERO)
                    .with_hydrophilic(true),
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
