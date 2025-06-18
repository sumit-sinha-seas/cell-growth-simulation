#pragma once

#include <array>

namespace cell_simulation {

class Particle {
public:
    Particle();
    ~Particle() = default;

    // Getters
    double getRadius() const { return radius_; }
    double getElasticModulus() const { return elastic_modulus_; }
    double getPoissonRatio() const { return poisson_ratio_; }
    double getReceptorConcentration() const { return receptor_concentration_; }
    double getLigandConcentration() const { return ligand_concentration_; }
    double getGrowthRate() const { return growth_rate_; }
    double getCellCycleTime() const { return cell_cycle_time_; }
    const std::array<double, 3>& getPosition() const { return position_; }
    const std::array<double, 3>& getVelocity() const { return velocity_; }
    const std::array<double, 3>& getForce() const { return force_; }

    // Setters
    void setRadius(double radius) { radius_ = radius; }
    void setElasticModulus(double modulus) { elastic_modulus_ = modulus; }
    void setPoissonRatio(double ratio) { poisson_ratio_ = ratio; }
    void setReceptorConcentration(double concentration) { receptor_concentration_ = concentration; }
    void setLigandConcentration(double concentration) { ligand_concentration_ = concentration; }
    void setGrowthRate(double rate) { growth_rate_ = rate; }
    void setCellCycleTime(double time) { cell_cycle_time_ = time; }
    void setPosition(const std::array<double, 3>& pos) { position_ = pos; }
    void setVelocity(const std::array<double, 3>& vel) { velocity_ = vel; }
    void setForce(const std::array<double, 3>& force) { force_ = force; }

    // Update methods
    void updatePosition(const std::array<double, 3>& delta_pos);
    void updateVelocity(const std::array<double, 3>& delta_vel);
    void updateForce(const std::array<double, 3>& delta_force);
    void grow(double delta_time);
    bool shouldDivide() const;

private:
    // Physical properties
    double radius_{0.0};
    double elastic_modulus_{0.0};
    double poisson_ratio_{0.0};
    double receptor_concentration_{0.0};
    double ligand_concentration_{0.0};
    double growth_rate_{0.0};
    double cell_cycle_time_{0.0};

    // Dynamic properties
    std::array<double, 3> position_{0.0, 0.0, 0.0};
    std::array<double, 3> velocity_{0.0, 0.0, 0.0};
    std::array<double, 3> force_{0.0, 0.0, 0.0};
};

} // namespace cell_simulation 