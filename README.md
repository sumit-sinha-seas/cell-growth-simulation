# 3D Cell Growth & Division Simulation

A sophisticated C++ simulation framework for modeling 3D cell growth, division, and mechanical interactions in biological tissues. This project provides a robust implementation of cell dynamics, including mechanical forces, growth mechanics, and cell division processes.

## 🌟 Features

### Core Functionality
- **3D Cell Dynamics**: Full 3D simulation of cell movement and interactions
- **Mechanical Interactions**:
  - Cell-cell repulsion forces
  - Adhesive forces based on receptor-ligand interactions
  - Pressure-dependent growth and division
- **Growth & Division**:
  - Realistic cell growth mechanics
  - Pressure-dependent division
  - Daughter cell inheritance of properties
- **Boundary Conditions**:
  - Periodic boundary conditions for infinite tissue simulation
  - Configurable simulation box dimensions

### Technical Features
- **Modern C++ Implementation**:
  - C++17 standard
  - Object-oriented design
  - Smart memory management
  - Exception handling
- **Build System**:
  - CMake-based build system
  - Cross-platform compatibility
  - Configurable compilation options
- **Data Output**:
  - LAMMPS-compatible trajectory files
  - Time-series data for analysis
  - Compatible with ParaView and OVITO visualization

## 🚀 Getting Started

### Prerequisites
- C++17 compatible compiler (GCC 7+, Clang 5+, or MSVC 2019+)
- CMake 3.10 or higher
- Git

### Installation

1. Clone the repository:
```bash
git clone https://github.com/sumit-sinha-seas/cell-growth-simulation.git
cd cell-growth-simulation
```

2. Create and enter build directory:
```bash
mkdir build && cd build
```

3. Configure and build:
```bash
cmake ..
make
```

### Running the Simulation

Basic usage:
```bash
./cell_simulation <initial_cells> <adhesion_strength> <simulation_time> <ensemble_size> <file_number>
```

Example:
```bash
./cell_simulation 200 0.1 10000 1 1
```

Parameters:
- `initial_cells`: Number of initial cells
- `adhesion_strength`: Cell-cell adhesion strength
- `simulation_time`: Total simulation time steps
- `ensemble_size`: Number of ensemble runs
- `file_number`: Output file identifier

## 📊 Output Files

The simulation generates several output files:
- `snapshot_data_*.dat`: Cell positions and properties at each snapshot
- `position_data_*.dat`: Time series of cell positions
- `cell_properties_*.dat`: Evolution of cell properties (elastic modulus, Poisson ratio, etc.)

## 🏗️ Project Structure

```
cell-growth-simulation/
├── CMakeLists.txt          # Build configuration
├── include/                # Header files
│   ├── cell_simulation.hpp # Main simulation class
│   ├── particle.hpp        # Cell particle class
│   └── simulation_parameters.hpp # Simulation parameters
├── src/                    # Source files
│   ├── main.cpp           # Entry point
│   ├── cell_simulation.cpp # Simulation implementation
│   └── particle.cpp       # Particle implementation
└── README.md              # This file
```

## 🔧 Configuration

Key simulation parameters can be adjusted in `include/simulation_parameters.hpp`:
- Box dimensions
- Physical parameters (viscosity, adhesion strength)
- Growth and division parameters
- Initial conditions

## 📈 Visualization

Output files can be visualized using:
- ParaView
- OVITO
- Custom Python scripts (examples provided in `scripts/`)

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Sumit Sinha** - *Initial work* - [sumit-sinha-seas](https://github.com/sumit-sinha-seas)

## 🙏 Acknowledgments

- Based on research in cell mechanics and tissue dynamics
- Inspired by various biological tissue simulation models
- Built with modern C++ best practices

## 📚 References

1. Cell Mechanics and Tissue Dynamics
2. Particle-based Simulation Methods
3. Biological Tissue Modeling

## 🔮 Future Work

- [ ] Add GPU acceleration
- [ ] Implement more complex cell behaviors
- [ ] Add Python bindings
- [ ] Improve visualization tools
- [ ] Add more comprehensive testing 