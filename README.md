# CLI_Fluid_Simulation
This is a repository for PostFix calculator using Stack data structure.
The source code was made for 5G00DM61-3003 Programming Languages 3.
## Concept
In this exercise, either different variations of the postfix calculator or different variations of the symbol balancing program are made. Specifications for variations are described below. Also below is a list of return requirements.

When doing the work, follow the order one-point job -> two-point job ->…, ie a job worth five points must contain all the required functions of levels 1 - 5. Thus, the work should be done incrementally, always step by step in this order. When moving to work worth the next point, don’t take any old qualities/code out of your work.

## Screenshot


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