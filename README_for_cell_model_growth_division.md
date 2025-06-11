# 3D Cell Growth & Division Simulation

This project implements a 3D particle-based simulation of cell growth, division, mechanical interaction, and death dynamics using C++. The model incorporates force laws, pressure-dependent growth, and stochastic death, simulating biologically inspired behavior in dense tissues or aggregates.

---

## 🧠 Features

- Cell-cell mechanical interactions: repulsion, adhesion, and pressure response
- Growth and division mechanics with receptor-ligand contact modeling
- Random cell death based on stress or crowding
- Periodic boundary conditions
- LAMMPS-style output compatible with ParaView or OVITO
- Time-logged `.dat` outputs for analysis

---

## 📦 File Description

- `cell_model_growth_division.cpp` — main simulation logic
- Output files:
  - `lv_Position_*.dat`: Position and radius of live cells
  - `dd_Position_*.dat`: Dead cells
  - `*_Cell_Properties_*.dat`: nu, E, c1, c2 parameters
  - `#Time_Record_*.dat`: Time taken to complete run

---

## 🛠️ Build Instructions (macOS/Linux)

### 1. Prerequisites

Ensure you have a C++ compiler like `g++` or `clang++`. Then compile:

```bash
g++ -O3 -std=c++17 cell_model_growth_division.cpp -o simulate
```

### 2. Run

```bash
./simulate <N0> <fad> <stopTime> <ensembleSize> <fileNum>
```

Example:

```bash
./simulate 200 0.1 10000 1 1
```

---

## 📊 Output

This simulation writes data to `.dat` files that can be loaded in scientific visualization tools like ParaView or OVITO. Each run records:

- Live and dead cell trajectories
- Cell property evolution
- Total simulation time

---

## 🧼 Notes

- This is a working prototype. A modular, object-oriented refactor is in progress.
- Tested on Linux and macOS with Clang/GCC.
- `<bits/stdc++.h>` should be removed for cross-platform compatibility.

---

## 📘 License

MIT License
