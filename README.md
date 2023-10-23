# rushydro

<p align="center">
    <img width="128" src="assets/icon-256.png" alt="rushydro icon">
</p>

Fluid simulation with particles

Try it on your browser! https://msakuta.github.io/rushydro/

## Gallery

![screenshot](images/screenshot00.png)


### Watermill

A demonstration of combined simulation of Lagrangean fluid and rigid bodies.

https://github.com/msakuta/rushydro/assets/2798715/97aecdee-ec10-4994-9bd3-7013f273e2cd


### Buoyancy

A solid object whose average density is lighter than the liquid will float on it.

https://github.com/msakuta/rushydro/assets/2798715/974de89f-6578-4add-983e-41e4f08c40ff


### Heat convection

Below is a heat convection simulation, where the bottom of the screen has a heat source that increases fluid temperature.
Higher tenperature particles exert more pressure. It is notable that convection (hot fluid floats, cool fluid sinks) naturally emerges.
Nearby particles can exchange heat, which realizes heat diffusion.
Note that the temperature is an attribute of each particle, so they will advect along with the particles.

https://github.com/msakuta/rushydro/assets/2798715/158b1270-1737-4e68-af25-41c63ea996f3





## What is this?

This is a sister project of [cfd-wasm](https://github.com/msakuta/cfd-wasm), a fluid simulator in Rust.
While cfd-wasm used _Eulerian method_ to discretize space, this project uses finite _Lagrangean method_, meaning we track attributes of particles.
It has more advantage that it can calculate varying density like gases or liquid surfaces.

This project has been inspired by [Sebastian Lague's video](https://youtu.be/rSKMYc1CQHE?si=4z0-JIuDQ7tOuDHR),
who did fluid simulation in Unity.
If you can do it in Unity, why not in Rust?


## How to run natively

```
cargo r
```


## How to build a WebAssembly version

Install trunk by 

```
cargo install --locked trunk
```

and run

```
trunk build --release
```

