#include "cell_simulation.hpp"
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    try {
        cell_simulation::CellGrowthSimulation simulation;
        
        // Initialize simulation
        simulation.initialize();
        
        // Run simulation
        simulation.run();
        
        // Cleanup
        simulation.cleanup();
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 