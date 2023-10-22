use eframe::epaint::{pos2, vec2, Vec2};
use rand::{thread_rng, Rng};
use std::{cell::Cell, collections::HashMap};

use super::{
    obstacle::{Environment, Obstacle},
    HashEntry, NeighborMode, NeighborPayload, Params, RusHydroApp, CELL_SIZE, MAX_HISTORY,
    PARTICLE_RADIUS, PARTICLE_RADIUS2,
};
use crate::measure_time;

pub(super) struct Particle {
    pub pos: Cell<Vec2>,
    pub velo: Cell<Vec2>,
    pub temp: Cell<f32>,
}

impl RusHydroApp {
    pub(super) fn reset(&mut self) {
        let mut rng = thread_rng();
        let particles = (0..self.num_particles)
            .map(|_| Particle {
                pos: Cell::new(vec2(
                    rng.gen_range(self.rect.min.x..self.rect.max.x),
                    rng.gen_range(self.rect.min.y..self.rect.max.y),
                )),
                velo: Cell::new(Vec2::ZERO),
                temp: Cell::new(0.),
            })
            .collect();
        self.particles = particles;
        // Force reinitialization of neighbor caches
        self.neighbor_payload = NeighborPayload::BruteForce;
    }

    fn repulsion_func(dist: f32, params: &Params, temp: f32) -> f32 {
        let radius = PARTICLE_RADIUS * (1.0 + 0.5 * temp);
        let f = (1. - dist / radius).powf(2.);
        1. / dist
            * (if f < params.surface_tension_threshold {
                (f - params.surface_tension_threshold) * params.surface_tension
            } else {
                f - params.surface_tension_threshold
            })
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
            let repulsion = delta
                * Self::repulsion_func(dist, params, particle_i.temp.get() + particle_j.temp.get());
            particle_i.velo.set(
                velo_i
                    + repulsion * params.repulsion_force
                    + (average_velo - velo_i) * params.viscosity,
            );
            particle_j.velo.set(
                velo_j - repulsion * params.repulsion_force
                    + (average_velo - velo_j) * params.viscosity,
            );
            let heat_conduction = 1. - dist / PARTICLE_RADIUS;
            let temp_i = particle_i.temp.get();
            let temp_j = particle_j.temp.get();
            let average_temp = (temp_i + temp_j) / 2.;
            particle_i
                .temp
                .set(temp_i + (average_temp - temp_i) * heat_conduction * 0.1);
            particle_j
                .temp
                .set(temp_j + (average_temp - temp_j) * heat_conduction * 0.1);
        }
    }

    fn update_speed_from_obstacles(
        obstacles: &mut [Obstacle],
        particle: &Particle,
        params: &Params,
    ) {
        let pos = particle.pos.get();
        let mut velo = particle.velo.get();
        for obstacle in obstacles {
            if !obstacle.is_hydrophilic() {
                continue;
            }
            if let Some(result) = obstacle.distance(pos) {
                if 0. < result.dist && result.dist < PARTICLE_RADIUS {
                    let impulse = result.normal
                        * Self::repulsion_func(result.dist, params, particle.temp.get());
                    velo = impulse;
                    obstacle.apply_impulse(-impulse, pos);
                }
            }
        }
        particle.velo.set(velo);
    }

    pub(super) fn update_particles(&mut self) {
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
                let (_, time) = measure_time(|| {
                    for (i, particle_i) in self.particles.iter().enumerate() {
                        for (j, particle_j) in self.particles.iter().enumerate() {
                            if i == j {
                                continue;
                            }
                            Self::update_speed(particle_i, particle_j, &self.params);
                        }
                        Self::update_speed_from_obstacles(
                            &mut self.obstacles,
                            particle_i,
                            &self.params,
                        );
                    }
                });
                self.time_history.push_back(time);
                if MAX_HISTORY < self.time_history.len() {
                    self.time_history.pop_front();
                }
            }
            NeighborPayload::HashMap(ref mut hash_map) => {
                let (_, hash_time) = measure_time(|| {
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

                let (_, update_time) = measure_time(|| {
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
                        Self::update_speed_from_obstacles(
                            &mut self.obstacles,
                            particle_i,
                            &self.params,
                        );
                    }
                });
                self.hash_time_history.push_back(super::HashMapTime {
                    hash: hash_time,
                    update: update_time,
                });
                if MAX_HISTORY < self.hash_time_history.len() {
                    self.hash_time_history.pop_front();
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
                    for (i, (particle_i, hash_entry)) in
                        self.particles.iter().zip(hash_table.iter_mut()).enumerate()
                    {
                        let pos = particle_i.pos.get();
                        let grid_pos = (
                            pos.x.div_euclid(CELL_SIZE) as i32,
                            pos.y.div_euclid(CELL_SIZE) as i32,
                        );
                        let cell_hash = hasher(grid_pos);
                        hash_entry.particle_idx = i;
                        hash_entry.cell_hash = cell_hash;
                    }
                });
                let (_, sort_time) = measure_time(|| {
                    hash_table.sort_unstable_by_key(|entry| entry.cell_hash);
                });
                let (_, find_time) = measure_time(|| {
                    for i in 0..self.particles.len() {
                        let first = hash_table.binary_search_by_key(&i, |entry| entry.cell_hash);
                        if let Ok(first) = first {
                            start_offsets[i] = 0;
                            for j in (0..=first).rev() {
                                if hash_table[j].cell_hash != i {
                                    start_offsets[i] = j + 1;
                                    break;
                                }
                            }
                        }
                    }
                });

                let (_, update_time) = measure_time(|| {
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
                        Self::update_speed_from_obstacles(
                            &mut self.obstacles,
                            particle_i,
                            &self.params,
                        );
                    }
                });
                self.sort_time_history.push_back(super::SortMapTime {
                    hash: hash_time,
                    find: find_time,
                    sort: sort_time,
                    update: update_time,
                });
                if MAX_HISTORY < self.sort_time_history.len() {
                    self.sort_time_history.pop_front();
                }
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

        for particle in &mut self.particles {
            let pos = particle.pos.get();
            let mut velo = particle.velo.get();
            velo.y -= self.gravity;
            let newpos = pos + velo;
            let mut croppos = pos2(
                newpos.x.min(self.rect.max.x).max(self.rect.min.x),
                newpos.y.min(self.rect.max.y).max(self.rect.min.y),
            );
            for obstacle in &mut self.obstacles {
                if let Some(result) = obstacle.hit_test(newpos) {
                    croppos -= result.normal * result.dist;
                    let impulse = result.normal
                        * result.normal.dot(velo - result.velo)
                        * (1. + self.restitution);
                    velo -= impulse;
                    obstacle.apply_impulse(impulse, newpos);
                }
            }
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
            if matches!(self.obstacle_select, Environment::Convection)
                && newpos.y < self.rect.min.y + 4.
                && newpos.x.abs() < self.rect.width() * 0.25
            {
                particle.temp.set((particle.temp.get() + 0.005).min(1.));
            } else {
                particle.temp.set((particle.temp.get() - 0.0005).max(0.));
            }
            match self.obstacle_select {
                Environment::Snake => {
                    if self.rect.max.x / 2. < croppos.x && croppos.y < self.rect.min.y + 1. {
                        croppos.x = -croppos.x;
                        croppos.y = self.rect.max.y - 1.;
                        velo.x = -velo.x;
                    }
                }
                Environment::Slope | Environment::WaterMill => {
                    if self.rect.max.x - PARTICLE_RADIUS < croppos.x {
                        croppos.x = self.rect.min.x + PARTICLE_RADIUS * 0.5;
                        croppos.y += self.rect.width() * Self::SLOPE_ANGLE.sin();
                        velo *= 0.5;
                    }
                }
                _ => {}
            }
            particle.pos.set(croppos.to_vec2());
            particle.velo.set(velo);
        }
    }
}
