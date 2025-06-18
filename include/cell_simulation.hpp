#pragma once

#include <vector>
#include <random>
#include <memory>
#include <string>
#include <fstream>
#include <stdexcept>

namespace cell_simulation {

// Forward declarations
class Particle;
class SimulationParameters;
class SimulationState;

class CellGrowthSimulation {
public:
    CellGrowthSimulation();
    ~CellGrowthSimulation() = default;

    // Prevent copying
    CellGrowthSimulation(const CellGrowthSimulation&) = delete;
    CellGrowthSimulation& operator=(const CellGrowthSimulation&) = delete;

    // Allow moving
    CellGrowthSimulation(CellGrowthSimulation&&) = default;
    CellGrowthSimulation& operator=(CellGrowthSimulation&&) = default;

    // Main simulation methods
    void initialize();
    void run();
    void cleanup();

private:
    // Core simulation components
    std::unique_ptr<SimulationParameters> params_;
    std::unique_ptr<SimulationState> state_;
    std::vector<std::unique_ptr<Particle>> particles_;

    // Random number generation
    std::random_device r_;
    std::seed_seq seed_;
    std::mt19937 gen_;

    // Internal methods
    void initializeRandomDistributions();
    void initializeParticles();
    void initializePositions();
    void calculateForces();
    void updatePositions();
    void handleParticleGrowth();
    void handleParticleDivision();
    void saveSnapshot() const;
    void savePositions() const;

    // Helper methods
    double applyPeriodicBoundaryConditions(double x) const;
    void calculatePairForces(size_t i, size_t j);
    double calculateRepulsiveForce(const Particle& p1, const Particle& p2) const;
    double calculateAdhesiveForce(const Particle& p1, const Particle& p2) const;
};

} // namespace cell_simulation 