# rushydro

Fluid simulation with particles

Try it on your browser! https://msakuta.github.io/rushydro/

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

