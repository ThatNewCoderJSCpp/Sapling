# Sapling

Physics simulations and sandboxes, built with [Fizmo](https://github.com/ThatNewCoderJSCpp/Fizmo)'s numerics, [Grove](https://github.com/ThatNewCoderJSCpp/Grove) for the window and [Walnut](https://github.com/ThatNewCoderJSCpp/Walnut) for its script. It's the physics companion to Arbor: Arbor graphs equations, Sapling runs worlds.

## Sections

| Section | Simulations so far |
|---|---|
| **Mechanics** | Double pendulum (adaptive Dormand–Prince, with a twin that starts a hair off to show chaos) |
| **Gravity** | Kepler orbit (Yoshida's symplectic integrator), orbit around a black hole (Schwarzschild geodesics, with the horizon, photon sphere, ISCO and the Newtonian orbit to compare) |
| **Quantum** | Tunneling through a barrier (Schrödinger's equation as a pair of real PDEs, Crank–Nicolson) |
| **Thermodynamics** | Heat in a rod (the heat equation with fixed end temperatures) |
| **Fluids and aero** | Flow past a spinning cylinder (streamlines, tracers and Kutta–Joukowski lift) |
| **Electromagnetism** | Field lines and equipotentials of moving point charges |

Each simulation has its settings as sliders, play, pause, step, reset and speed controls, and live readouts.

## Building

Sapling needs CMake 3.18 or newer, a C++17 compiler, Ninja (recommended) and the same system libraries as Grove (X11, Xft, fontconfig, FreeType, HarfBuzz, OpenSSL, zlib, libjpeg and the Vulkan headers on Linux).

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
build/bin/sapling
```

CMake fetches Fizmo (`SAPLING_FIZMO_TAG`), Grove (`SAPLING_GROVE_TAG`) and Walnut (`SAPLING_WALNUT_TAG`). Sapling declares Fizmo itself, before Grove, so the Fizmo it asks for is the one Walnut and Grove build with too. It needs Fizmo's `numerics` component (v0.2.38 or newer). To build against a local Fizmo checkout instead, pass `-DFETCHCONTENT_SOURCE_DIR_FIZMO=/path/to/Fizmo`.

| Option | Default | What it does |
|---|---|---|
| `SAPLING_BUILD_TESTS` | ON | builds `sapling_tests` |
| `SAPLING_BUNDLE` | ON | runs `sapling-host --bundle` to make `sapling` with its window and compiled script built in |
| `SAPLING_PACKAGE` | OFF | makes release packages with CPack |
| `SAPLING_WITH_WALNUT` | ON | builds Walnut alongside Grove for the script |

`sapling-host` runs the UI straight from the `UI` folder, which is handy while working on it; set `SAPLING_UI_DIR` to point it somewhere else.

## How it's put together

- `Include/sapling`, `Source`: the core library. `Simulation` is the base class (reset, step, draw, readouts and parameters); each section's simulations are in `Source/Simulations`. `Canvas` is a small anti-aliased software rasterizer the host draws each frame with.
- `Tools/sapling.cpp`: the host program. It owns the running simulation, steps it on a timer and hands the pixels to Grove's `<canvas>`.
- `UI`: the window (`sapling.grove`), its styles (`sapling.bark`) and its script (`sapling.wal`).
- `Tests/sapling_tests.cpp`: every simulation runs and draws, plus physics checks (energy conservation, a one-year circular orbit, periapsis precession and plunging near a black hole, conserved probability in tunneling, the rod's steady state, Kutta–Joukowski lift).

Adding a simulation: write a class deriving from `Simulation` in `Source/Simulations`, add a `make_...` function to `simulations.hpp`, list it in `simulations()` and `make_simulation()` in `Source/simulation.cpp`, and add a test.

## Branches and releases

Work happens on `dev`, and `main` is fast-forwarded to it. Pushes to either branch run CI (`.github/workflows/ci.yml`). Pushing a `v*` tag runs `.github/workflows/release.yml`, which builds Linux (x86_64, arm64) and Windows (x86_64, i686, arm64 with llvm-mingw) packages and attaches them to the release. The tag has to match `VERSION` in `CMakeLists.txt`.
