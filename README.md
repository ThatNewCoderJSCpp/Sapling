# Sapling

A physics and chemistry sandbox, built with [Fizmo](https://github.com/ThatNewCoderJSCpp/Fizmo)'s numerics and chemistry, [Grove](https://github.com/ThatNewCoderJSCpp/Grove) for the window and [Walnut](https://github.com/ThatNewCoderJSCpp/Walnut) for its script. It's the physics companion to Arbor: Arbor graphs equations, Sapling runs worlds.

The window has two sections, picked at the top.

**Physics** is being rebuilt. For now it has the stage, play, pause, step, reset and speed controls, and the settings and readouts panels, ready for whatever is running in it.

**Chemistry** has these tabs:
- **Bonding sandbox**: covalent, ionic and metallic bonding, each in 2D and 3D.
  - *Covalent*: pick nonmetals and metalloids (aluminum and beryllium too) from a periodic table and build molecules yourself. Click to place atoms, drag between atoms to bond them, and click a bond to make it double or triple. Fill H, Tidy (a VSEPR shape in 3D) and Fix Lewis help out, and polyatomic ions and common molecules load as templates. The side panel tracks the name, stability verdict, Lewis structure (lone pairs, formal charges, broken octets, radicals), resonance and average bond orders, the shape around each atom, polarity and dipole, and ΔH. Click a bond to see its order (whole or fractional), average length, bond energy, ΔH, ΔS and ΔG to break it at the chosen temperature and pressure, how many are broken, its polarity and its stretching frequency. Lone pairs, formal and partial charges, bond dipoles and the molecular dipole can be shown on the drawing.
  - *Ionic*: pick a metal and a nonmetal (or type any ions, polyatomic ones included). The crystal is built in its likely structure and assembles itself from scattered ions. The panel gives the lattice energy four ways with the Born–Haber cycle, radii, coordination, density, solubility and more.
  - *Metallic*: pick up to four metals and their proportions, or a preset (brass, bronze, steel, solder, gold). The panel names the alloy and checks the Hume-Rothery rules; the drawing shows the lattice with its electron sea, which drifts when an electric field is turned on.
  - Temperature and pressure (any units) change bond breaking, and melt or boil lattices past their melting and boiling points.
  - **Save** writes the sandbox to `.sapling/chemistry/covalent`, `ionic` or `metallic` in the folder Sapling runs from (`SAPLING_SAVE_DIR` moves it; on Android it lives in the app's data folder). `.sapling/` is in `.gitignore`.
- **Ions & names**: 157 polyatomic ions you can search and filter, each with a 3D or Lewis structure, its other names and its acid. A namer turns formulas into names and names into formulas (ionic compounds, acids, hydrates, molecular compounds, organics, carbohydrates), and two ions combine into a compound.
- **Reactions**: predict the products from only the reactants (combustion, synthesis, decomposition, single and double replacement, neutralization, equilibria with ⇌), balance an equation with your choice of arrow, see ionic equations, redox changes and ΔH, ΔS, ΔG and K. A stoichiometry solver takes masses, moles, gas volumes or solution volumes of any species and gives the limiting reactant and every amount.
- **Solvers**: over 60 equations (gases, solutions, acids and bases, thermochemistry, kinetics, electrochemistry, light). Fill in what you know in any units and solve for the rest. Also a unit converter, percent composition, and empirical and molecular formulas.
- **Phase diagrams**: from a substance name, formula, SMILES, element or anything saved in the bonding sandbox. You get the 2D pressure–temperature diagram with a draggable state point, and the 3D pressure–specific volume–temperature surface. Every property can be edited to make your own diagram.
- **Radial probability**: P(r) for each subshell of any element or ion, with the total, cumulative curves, peaks and nodes.
- **Electron clouds**: any orbital or whole atom as a 3D point cloud or a 2D density slice, colored by the sign of ψ.
- **Elements**: a periodic table with each element's shells, orbital-box diagram and properties.
- **Units**: the units every result is shown in. Inputs with units have their own dropdowns.

## Building

Sapling needs CMake 3.18 or newer, a C++17 compiler, Ninja (recommended) and the same system libraries as Grove (X11, Xft, fontconfig, FreeType, HarfBuzz, OpenSSL, zlib, libjpeg and the Vulkan headers on Linux).

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
build/bin/sapling
```

CMake fetches Fizmo (`SAPLING_FIZMO_TAG`), Grove (`SAPLING_GROVE_TAG`) and Walnut (`SAPLING_WALNUT_TAG`). Sapling declares Fizmo itself, before Grove, so the Fizmo it asks for is the one Walnut and Grove build with too. It needs Fizmo's `numerics` and `chemistry` components (v0.3.5 or newer). To build against a local Fizmo checkout instead, pass `-DFETCHCONTENT_SOURCE_DIR_FIZMO=/path/to/Fizmo`.

| Option | Default | What it does |
|---|---|---|
| `SAPLING_BUILD_TESTS` | ON | builds `sapling_tests` |
| `SAPLING_BUNDLE` | ON | runs `sapling-host --bundle` to make `sapling` with its window and compiled script built in |
| `SAPLING_PACKAGE` | OFF | makes release packages with CPack |
| `SAPLING_WITH_WALNUT` | ON | builds Walnut alongside Grove for the script |

`sapling-host` runs the UI straight from the `UI` folder, which is handy while working on it; set `SAPLING_UI_DIR` to point it somewhere else.

## How it's put together

- `Include/sapling`, `Source`: the core library. `Simulation` is the base class for anything the sandbox runs (reset, step, draw, readouts and parameters). `Canvas` is a small anti-aliased software rasterizer the host draws each frame with: lines, shapes, polygons, shaded spheres and a stroke font with Greek letters, sub- and superscripts.
- `Include/sapling/chemistry.hpp`, `Source/chemistry`: the chemistry views (`ChemView`) the host draws into canvases. There's the bonding sandbox, structure, phase diagram, radial probability, electron cloud and element views, plus a small 3D renderer, plot axes, unit preferences, the JSON save files and `chemistry_tool` for the naming, balancing, prediction and solver panels. All the chemistry itself comes from Fizmo's `chemistry` component.
- `Tools/sapling.cpp`: the host program. It owns what's running, steps it on a timer and hands the pixels to Grove's `<canvas>`; with nothing running it draws an empty grid. It also routes mouse and keys to the chemistry canvases, and gives the scripts the `sapling.chem.*` functions.
- `UI`: the window (`sapling.grove`), its styles (`sapling.bark`, `chemistry.bark`), its scripts (`sapling.wal`, `chemistry.wal`) and one file per chemistry tab in `UI/chemistry`.
- `Tests/sapling_tests.cpp`: the canvas, the `Simulation` base class (stepping, parameters and their ranges, readouts, drawing), number formatting, and the chemistry views, tools, JSON and saving.

## Branches and releases

Work happens on `dev`, and `main` is fast-forwarded to it. Pushes to either branch run CI (`.github/workflows/ci.yml`). Pushing a `v*` tag runs `.github/workflows/release.yml`, which builds Linux (x86_64, arm64), macOS (arm64, x86_64), Android (one APK) and Windows (x86_64, i686, arm64 with llvm-mingw) packages and attaches them to the release. The tag has to match `VERSION` in `CMakeLists.txt`.

### macOS

`Sapling-<version>-macos-arm64.tar.gz` is for Apple Silicon Macs and `macos-x86_64` for Intel Macs (macOS 11 or newer). Everything Sapling needs that macOS doesn't come with (FreeType, HarfBuzz, fontconfig, libjpeg and what they use) is in the package's `lib` folder, so nothing else has to be installed. With `-DSAPLING_PACKAGE=ON` on macOS, the build copies those libraries out of Homebrew into `build/bin/lib` (`cmake/macos_bundle_libs.cmake`) and points the programs at them.

The download isn't signed by Apple, so macOS refuses to open it at first. Remove the download flag once, then run it:

```sh
xattr -dr com.apple.quarantine Sapling-*
./Sapling-*/sapling
```

### Android

`Sapling-<version>-android.apk` runs on phones, tablets and Chromebooks. Copy it to the phone and open it there (Android asks once to allow installing apps from that source), or install it over USB with `adb install -r Sapling-*-android.apk`. The release workflow's `.github/scripts/android-apk.sh` builds it: a desktop `grove` and `walnut`, then `libsapling_android.so` and Grove's Android kit for each ABI with the Android NDK, then `grove build --target android` compiles the script for the phone and makes the APK. Signing uses the `SAPLING_KEYSTORE_BASE64`, `SAPLING_KEYSTORE_PASS` and `SAPLING_KEY_ALIAS` secrets when they're set; without them each build gets a new debug key, so an older install has to be uninstalled first.

