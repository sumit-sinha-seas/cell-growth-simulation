#include "particle.hpp"
#include <cmath>

namespace cell_simulation {

Particle::Particle() = default;

void Particle::updatePosition(const std::array<double, 3>& delta_pos) {
    for (size_t i = 0; i < 3; ++i) {
        position_[i] += delta_pos[i];
    }
}

void Particle::updateVelocity(const std::array<double, 3>& delta_vel) {
    for (size_t i = 0; i < 3; ++i) {
        velocity_[i] += delta_vel[i];
    }
}

void Particle::updateForce(const std::array<double, 3>& delta_force) {
    for (size_t i = 0; i < 3; ++i) {
        force_[i] += delta_force[i];
    }
}

void Particle::grow(double delta_time) {
    radius_ += growth_rate_ * delta_time;
    cell_cycle_time_ += delta_time;
}

bool Particle::shouldDivide() const {
    return radius_ >= 2.0 * mitotic_radius_;
}

} // namespace cell_simulation 