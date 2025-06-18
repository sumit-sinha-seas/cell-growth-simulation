#pragma once

namespace cell_simulation {

struct SimulationParameters {
    // Simulation dimensions
    static constexpr int SPATIAL_DIMENSIONS = 3;
    
    // Time parameters
    int start_time{0};
    int stop_time{0};
    int time_step{10};
    int snapshot_interval{10};
    
    // Box parameters
    double box_dimension{0.0};
    double packing_fraction{0.0};
    
    // Physical parameters
    double cell_adhesion_strength{0.0};
    double ecm_viscosity{5e-3};
    double mitotic_radius{5.0};
    double critical_pressure{1e-4};
    
    // Distribution parameters
    double radius_mean{4.5};
    double radius_stddev{0.5};
    double poisson_ratio_mean{0.5};
    double poisson_ratio_stddev{0.02};
    double elastic_modulus_mean{1e-3};
    double elastic_modulus_stddev{1e-4};
    double receptor_concentration_mean{0.9};
    double receptor_concentration_stddev{0.02};
    
    // Growth parameters
    double growth_rate_stddev{1e-5};
    double cell_cycle_time{54000};
    
    // Initial conditions
    int initial_particle_count{0};
    int min_data_points{170};
    
    // Derived parameters
    double growth_rate_mean() const {
        return (2 * M_PI * mitotic_radius * mitotic_radius * mitotic_radius) / (3 * cell_cycle_time);
    }
    
    double time_step_size() const {
        return time_step * 1e-6;
    }
};

} // namespace cell_simulation 