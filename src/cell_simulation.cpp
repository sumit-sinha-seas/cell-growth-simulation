#include "cell_simulation.hpp"
#include "particle.hpp"
#include "simulation_parameters.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

namespace cell_simulation {

CellGrowthSimulation::CellGrowthSimulation()
    : seed_{r_(), r_(), r_(), r_(), r_()}
    , gen_(seed_)
    , params_(std::make_unique<SimulationParameters>())
    , state_(std::make_unique<SimulationState>()) {
    initializeRandomDistributions();
}

void CellGrowthSimulation::initialize() {
    initializeParticles();
    initializePositions();
}

void CellGrowthSimulation::run() {
    for (int t = params_->start_time; t <= params_->stop_time; t += params_->time_step) {
        calculateForces();
        updatePositions();
        handleParticleGrowth();
        handleParticleDivision();
        
        if (t % params_->snapshot_interval == 0) {
            saveSnapshot();
            savePositions();
        }
    }
}

void CellGrowthSimulation::cleanup() {
    particles_.clear();
}

void CellGrowthSimulation::initializeRandomDistributions() {
    // Initialize distributions with parameters from params_
    radius_dist_ = std::normal_distribution<double>(params_->radius_mean, params_->radius_stddev);
    poisson_ratio_dist_ = std::normal_distribution<double>(params_->poisson_ratio_mean, params_->poisson_ratio_stddev);
    elastic_modulus_dist_ = std::normal_distribution<double>(params_->elastic_modulus_mean, params_->elastic_modulus_stddev);
    concentration_dist_ = std::normal_distribution<double>(params_->receptor_concentration_mean, params_->receptor_concentration_stddev);
    growth_rate_dist_ = std::normal_distribution<double>(params_->growth_rate_mean(), params_->growth_rate_stddev);
}

void CellGrowthSimulation::initializeParticles() {
    particles_.reserve(params_->initial_particle_count);
    
    for (int i = 0; i < params_->initial_particle_count; ++i) {
        auto particle = std::make_unique<Particle>();
        
        // Set physical properties
        particle->setRadius(radius_dist_(gen_));
        particle->setElasticModulus(elastic_modulus_dist_(gen_));
        particle->setPoissonRatio(poisson_ratio_dist_(gen_));
        particle->setReceptorConcentration(concentration_dist_(gen_));
        particle->setLigandConcentration(concentration_dist_(gen_));
        particle->setGrowthRate(growth_rate_dist_(gen_));
        
        particles_.push_back(std::move(particle));
    }
}

void CellGrowthSimulation::initializePositions() {
    std::uniform_real_distribution<double> pos_dist(-params_->box_dimension/2.0, params_->box_dimension/2.0);
    
    for (auto& particle : particles_) {
        std::array<double, 3> position;
        for (int i = 0; i < 3; ++i) {
            position[i] = pos_dist(gen_);
        }
        particle->setPosition(position);
    }
}

void CellGrowthSimulation::calculateForces() {
    // Reset all forces to zero
    for (auto& particle : particles_) {
        particle->setForce({0.0, 0.0, 0.0});
    }
    
    // Calculate forces between all pairs of particles
    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            calculatePairForces(i, j);
        }
    }
}

void CellGrowthSimulation::calculatePairForces(size_t i, size_t j) {
    const auto& p1 = particles_[i];
    const auto& p2 = particles_[j];
    
    // Calculate distance vector with periodic boundary conditions
    std::array<double, 3> r;
    for (int k = 0; k < 3; ++k) {
        r[k] = applyPeriodicBoundaryConditions(p1->getPosition()[k] - p2->getPosition()[k]);
    }
    
    double distance_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    double distance = std::sqrt(distance_squared);
    double sum_radii = p1->getRadius() + p2->getRadius();
    
    if (distance < sum_radii) {
        double overlap = sum_radii - distance;
        
        // Calculate forces
        double F_rep = calculateRepulsiveForce(*p1, *p2);
        double F_ad = calculateAdhesiveForce(*p1, *p2);
        double F_total = F_rep + F_ad;
        
        // Update forces
        std::array<double, 3> force;
        for (int k = 0; k < 3; ++k) {
            force[k] = F_total * r[k] / distance;
        }
        
        p1->updateForce(force);
        for (int k = 0; k < 3; ++k) {
            force[k] = -force[k];
        }
        p2->updateForce(force);
    }
}

double CellGrowthSimulation::calculateRepulsiveForce(const Particle& p1, const Particle& p2) const {
    double E_eff = 2.0 * p1.getElasticModulus() * p2.getElasticModulus() / 
                  (p1.getElasticModulus() + p2.getElasticModulus());
    double nu_eff = 2.0 * p1.getPoissonRatio() * p2.getPoissonRatio() / 
                   (p1.getPoissonRatio() + p2.getPoissonRatio());
    
    return (4.0/3.0) * E_eff * std::sqrt(p1.getRadius() * p2.getRadius() / 
           (p1.getRadius() + p2.getRadius())) * 
           std::pow(overlap_, 1.5) * (1.0 - nu_eff);
}

double CellGrowthSimulation::calculateAdhesiveForce(const Particle& p1, const Particle& p2) const {
    double c_factor = p1.getReceptorConcentration() * p2.getLigandConcentration() + 
                     p2.getReceptorConcentration() * p1.getLigandConcentration();
    return params_->cell_adhesion_strength * c_factor * overlap_;
}

void CellGrowthSimulation::updatePositions() {
    double dt = params_->time_step_size();
    
    for (auto& particle : particles_) {
        // Update velocities
        std::array<double, 3> delta_vel;
        for (int i = 0; i < 3; ++i) {
            delta_vel[i] = dt * particle->getForce()[i];
        }
        particle->updateVelocity(delta_vel);
        
        // Update positions
        std::array<double, 3> delta_pos;
        for (int i = 0; i < 3; ++i) {
            delta_pos[i] = dt * particle->getVelocity()[i];
        }
        particle->updatePosition(delta_pos);
        
        // Apply periodic boundary conditions
        std::array<double, 3> pos = particle->getPosition();
        for (int i = 0; i < 3; ++i) {
            pos[i] = applyPeriodicBoundaryConditions(pos[i]);
        }
        particle->setPosition(pos);
    }
}

void CellGrowthSimulation::handleParticleGrowth() {
    double dt = params_->time_step_size();
    
    for (auto& particle : particles_) {
        particle->grow(dt);
    }
}

void CellGrowthSimulation::handleParticleDivision() {
    std::vector<std::unique_ptr<Particle>> new_particles;
    
    for (auto& particle : particles_) {
        if (particle->shouldDivide()) {
            // Create daughter particle
            auto daughter = std::make_unique<Particle>(*particle);
            
            // Adjust properties for division
            particle->setRadius(particle->getRadius() / 2.0);
            daughter->setRadius(particle->getRadius());
            
            // Reset cell cycle time
            particle->setCellCycleTime(0.0);
            daughter->setCellCycleTime(0.0);
            
            // Add to new particles
            new_particles.push_back(std::move(daughter));
        }
    }
    
    // Add new particles to the simulation
    particles_.insert(particles_.end(), 
                     std::make_move_iterator(new_particles.begin()),
                     std::make_move_iterator(new_particles.end()));
}

void CellGrowthSimulation::saveSnapshot() const {
    std::string filename = "snapshot_data_" + std::to_string(params_->start_time) + ".dat";
    std::ofstream outfile(filename);
    
    if (!outfile) {
        throw std::runtime_error("Failed to open snapshot data file for writing");
    }
    
    outfile << std::fixed << std::setprecision(6);
    for (const auto& particle : particles_) {
        const auto& pos = particle->getPosition();
        outfile << std::setw(10) << pos[0] << " "
                << std::setw(10) << pos[1] << " "
                << std::setw(10) << pos[2] << " "
                << std::setw(10) << particle->getRadius() << std::endl;
    }
}

void CellGrowthSimulation::savePositions() const {
    std::string filename = "position_data_" + std::to_string(params_->start_time) + ".dat";
    std::ofstream outfile(filename);
    
    if (!outfile) {
        throw std::runtime_error("Failed to open position data file for writing");
    }
    
    outfile << std::fixed << std::setprecision(6);
    for (const auto& particle : particles_) {
        const auto& pos = particle->getPosition();
        outfile << std::setw(10) << pos[0] << " "
                << std::setw(10) << pos[1] << " "
                << std::setw(10) << pos[2] << std::endl;
    }
}

double CellGrowthSimulation::applyPeriodicBoundaryConditions(double x) const {
    if (x > params_->box_dimension / 2.0) {
        return x - params_->box_dimension;
    } else if (x < -params_->box_dimension / 2.0) {
        return x + params_->box_dimension;
    }
    return x;
}

} // namespace cell_simulation 