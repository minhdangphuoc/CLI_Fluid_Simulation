# CLI_Fluid_Simulation
## Concept
- Basic concept of Fluid Simulation.
- Run on terminal.
- Textbase UI, Enhanced ASCII. -> Run on linux only for safety. Linux terminal supports full of symbols.
## Running
![Running Gif](./img/running.gif)

## Manual
This project using meson to build, so you have to download and install it to your computer. Meson Build

## Build
meson setup builddir //  Initialize the build
cd builddir // Move to build dirrectory
ninja && ./main // Run built
Change name of built
Change this code in meson.build for changing the name.

// meson.build
...
executable('main', 'main.cpp', link_with : lib) // change main to another name
...
## Issues
- Not found yet.
## References
[Fluid simulation for dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html)