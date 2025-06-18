#include <bits/stdc++.h>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <cmath>

// Remove using namespace std as it's considered bad practice
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::setprecision;
using std::fixed;
using std::setw;
using std::setfill;
using std::left;
using std::right;
using std::to_string;
using std::runtime_error;
using std::unique_ptr;
using std::make_unique;

// Constants namespace to avoid global variables
namespace Constants {
    constexpr int SP_DIM = 3;
    constexpr double EPS = 1e-30;
    constexpr double SIGMA = 0.01;
    constexpr double GAMMA_MAX = 1e-4;
    constexpr double RM = 5.0;
    constexpr double PC = 1e-4;
    constexpr double GR_RATE_RD_STDV = 1e-5;
}

// Main simulation class
class CellGrowthSimulation {
private:
    // Random number generation
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r()};
    std::mt19937 gen;

    // Simulation parameters
    struct Parameters {
        int N0{0};                    // Initial number of particles
        double boxDimension{0.0};     // Dimension of cubic box
        double fad{0};                // Cell-cell adhesion strength
        int startTime{0};
        int stopTime{0};
        int timeStep{10};
        int SC_startTime{10};
        int dynamicVectorDim{0};
        int id_lv{1}, id_DV{2}, id_dd{3};
        int minDataPoints{170};
        double eta{5e-3};            // ECM viscosity
        double l_mean{0};
        double l_stdv{25};
        double R_mean{4.5};
        double R_stdv{0.5};
        double mu_mean{0.045};       // Self propulsion mobility
        double nu_mean{0.5};
        double nu_stdv{0.02};        // Poisson ratio
        double E_mean{1e-3};
        double E_stdv{1e-4};         // Elastic moduli
        double c_mean{0.9};
        double c_stdv{0.02};
        double CC_time{54000};
        double rv{(2 * M_PI * std::pow(RM, 3)) / (3 * CC_time)};
        double kdt{timeStep * 1e-6};
        double scaleFactor{0.0};
        double flag_1{-50.0};
        double flag_2{-10.0};
    } params;

    // State variables
    struct State {
        double phi{0};               // Packing fraction
        int numParticles{0};         // Total number of cells
        int dd_numParticles{0};
        int DV_numParticles{0};
        double dx{0}, dy{0}, dz{0}, DIST{0}, DIST_sq{0};
        double R_cutoff{0.0};
        double Rd{0.0};
        double delta_Rd{0.0};
        double Rdtr{0.0};
        double PRESSURE{0.0};
        double E{0.0}, nu{0.0}, c1{0.0}, c2{0.0};
        double max_DIA{0.0};
        double overlap{0.0};
    } state;

    // Data containers
    struct DataContainers {
        vector<double> radius_list;
        vector<double> pos_x;
        vector<double> pos_y;
        vector<double> pos_z;
        vector<double> vel_x;
        vector<double> vel_y;
        vector<double> vel_z;
        vector<double> force_x;
        vector<double> force_y;
        vector<double> force_z;
        vector<double> elastic_moduli;
        vector<double> poisson_ratio;
        vector<double> receptor_concentration;
        vector<double> ligand_concentration;
        vector<double> growth_rate;
        vector<double> growth_rate_daughter;
        vector<double> cell_cycle_time;
        vector<double> cell_cycle_time_daughter;
        vector<double> cell_cycle_time_initial;
        vector<double> cell_cycle_time_initial_daughter;
        vector<double> cell_cycle_time_initial_1;
        vector<double> cell_cycle_time_initial_1_daughter;
        vector<double> cell_cycle_time_initial_2;
        vector<double> cell_cycle_time_initial_2_daughter;
        vector<double> cell_cycle_time_initial_3;
        vector<double> cell_cycle_time_initial_3_daughter;
        vector<double> cell_cycle_time_initial_4;
        vector<double> cell_cycle_time_initial_4_daughter;
        vector<double> cell_cycle_time_initial_5;
        vector<double> cell_cycle_time_initial_5_daughter;
        vector<double> cell_cycle_time_initial_6;
        vector<double> cell_cycle_time_initial_6_daughter;
        vector<double> cell_cycle_time_initial_7;
        vector<double> cell_cycle_time_initial_7_daughter;
        vector<double> cell_cycle_time_initial_8;
        vector<double> cell_cycle_time_initial_8_daughter;
        vector<double> cell_cycle_time_initial_9;
        vector<double> cell_cycle_time_initial_9_daughter;
        vector<double> cell_cycle_time_initial_10;
        vector<double> cell_cycle_time_initial_10_daughter;
        vector<double> cell_cycle_time_initial_11;
        vector<double> cell_cycle_time_initial_11_daughter;
        vector<double> cell_cycle_time_initial_12;
        vector<double> cell_cycle_time_initial_12_daughter;
        vector<double> cell_cycle_time_initial_13;
        vector<double> cell_cycle_time_initial_13_daughter;
        vector<double> cell_cycle_time_initial_14;
        vector<double> cell_cycle_time_initial_14_daughter;
        vector<double> cell_cycle_time_initial_15;
        vector<double> cell_cycle_time_initial_15_daughter;
        vector<double> cell_cycle_time_initial_16;
        vector<double> cell_cycle_time_initial_16_daughter;
        vector<double> cell_cycle_time_initial_17;
        vector<double> cell_cycle_time_initial_17_daughter;
        vector<double> cell_cycle_time_initial_18;
        vector<double> cell_cycle_time_initial_18_daughter;
        vector<double> cell_cycle_time_initial_19;
        vector<double> cell_cycle_time_initial_19_daughter;
        vector<double> cell_cycle_time_initial_20;
        vector<double> cell_cycle_time_initial_20_daughter;
        vector<double> cell_cycle_time_initial_21;
        vector<double> cell_cycle_time_initial_21_daughter;
        vector<double> cell_cycle_time_initial_22;
        vector<double> cell_cycle_time_initial_22_daughter;
        vector<double> cell_cycle_time_initial_23;
        vector<double> cell_cycle_time_initial_23_daughter;
        vector<double> cell_cycle_time_initial_24;
        vector<double> cell_cycle_time_initial_24_daughter;
        vector<double> cell_cycle_time_initial_25;
        vector<double> cell_cycle_time_initial_25_daughter;
        vector<double> cell_cycle_time_initial_26;
        vector<double> cell_cycle_time_initial_26_daughter;
        vector<double> cell_cycle_time_initial_27;
        vector<double> cell_cycle_time_initial_27_daughter;
        vector<double> cell_cycle_time_initial_28;
        vector<double> cell_cycle_time_initial_28_daughter;
        vector<double> cell_cycle_time_initial_29;
        vector<double> cell_cycle_time_initial_29_daughter;
        vector<double> cell_cycle_time_initial_30;
        vector<double> cell_cycle_time_initial_30_daughter;
        vector<double> cell_cycle_time_initial_31;
        vector<double> cell_cycle_time_initial_31_daughter;
        vector<double> cell_cycle_time_initial_32;
        vector<double> cell_cycle_time_initial_32_daughter;
        vector<double> cell_cycle_time_initial_33;
        vector<double> cell_cycle_time_initial_33_daughter;
        vector<double> cell_cycle_time_initial_34;
        vector<double> cell_cycle_time_initial_34_daughter;
        vector<double> cell_cycle_time_initial_35;
        vector<double> cell_cycle_time_initial_35_daughter;
        vector<double> cell_cycle_time_initial_36;
        vector<double> cell_cycle_time_initial_36_daughter;
        vector<double> cell_cycle_time_initial_37;
        vector<double> cell_cycle_time_initial_37_daughter;
        vector<double> cell_cycle_time_initial_38;
        vector<double> cell_cycle_time_initial_38_daughter;
        vector<double> cell_cycle_time_initial_39;
        vector<double> cell_cycle_time_initial_39_daughter;
        vector<double> cell_cycle_time_initial_40;
        vector<double> cell_cycle_time_initial_40_daughter;
        vector<double> cell_cycle_time_initial_41;
        vector<double> cell_cycle_time_initial_41_daughter;
        vector<double> cell_cycle_time_initial_42;
        vector<double> cell_cycle_time_initial_42_daughter;
        vector<double> cell_cycle_time_initial_43;
        vector<double> cell_cycle_time_initial_43_daughter;
        vector<double> cell_cycle_time_initial_44;
        vector<double> cell_cycle_time_initial_44_daughter;
        vector<double> cell_cycle_time_initial_45;
        vector<double> cell_cycle_time_initial_45_daughter;
        vector<double> cell_cycle_time_initial_46;
        vector<double> cell_cycle_time_initial_46_daughter;
        vector<double> cell_cycle_time_initial_47;
        vector<double> cell_cycle_time_initial_47_daughter;
        vector<double> cell_cycle_time_initial_48;
        vector<double> cell_cycle_time_initial_48_daughter;
        vector<double> cell_cycle_time_initial_49;
        vector<double> cell_cycle_time_initial_49_daughter;
        vector<double> cell_cycle_time_initial_50;
        vector<double> cell_cycle_time_initial_50_daughter;
        vector<double> cell_cycle_time_initial_51;
        vector<double> cell_cycle_time_initial_51_daughter;
        vector<double> cell_cycle_time_initial_52;
        vector<double> cell_cycle_time_initial_52_daughter;
        vector<double> cell_cycle_time_initial_53;
        vector<double> cell_cycle_time_initial_53_daughter;
        vector<double> cell_cycle_time_initial_54;
        vector<double> cell_cycle_time_initial_54_daughter;
        vector<double> cell_cycle_time_initial_55;
        vector<double> cell_cycle_time_initial_55_daughter;
        vector<double> cell_cycle_time_initial_56;
        vector<double> cell_cycle_time_initial_56_daughter;
        vector<double> cell_cycle_time_initial_57;
        vector<double> cell_cycle_time_initial_57_daughter;
        vector<double> cell_cycle_time_initial_58;
        vector<double> cell_cycle_time_initial_58_daughter;
        vector<double> cell_cycle_time_initial_59;
        vector<double> cell_cycle_time_initial_59_daughter;
        vector<double> cell_cycle_time_initial_60;
        vector<double> cell_cycle_time_initial_60_daughter;
        vector<double> cell_cycle_time_initial_61;
        vector<double> cell_cycle_time_initial_61_daughter;
        vector<double> cell_cycle_time_initial_62;
        vector<double> cell_cycle_time_initial_62_daughter;
        vector<double> cell_cycle_time_initial_63;
        vector<double> cell_cycle_time_initial_63_daughter;
        vector<double> cell_cycle_time_initial_64;
        vector<double> cell_cycle_time_initial_64_daughter;
        vector<double> cell_cycle_time_initial_65;
        vector<double> cell_cycle_time_initial_65_daughter;
        vector<double> cell_cycle_time_initial_66;
        vector<double> cell_cycle_time_initial_66_daughter;
        vector<double> cell_cycle_time_initial_67;
        vector<double> cell_cycle_time_initial_67_daughter;
        vector<double> cell_cycle_time_initial_68;
        vector<double> cell_cycle_time_initial_68_daughter;
        vector<double> cell_cycle_time_initial_69;
        vector<double> cell_cycle_time_initial_69_daughter;
        vector<double> cell_cycle_time_initial_70;
        vector<double> cell_cycle_time_initial_70_daughter;
        vector<double> cell_cycle_time_initial_71;
        vector<double> cell_cycle_time_initial_71_daughter;
        vector<double> cell_cycle_time_initial_72;
        vector<double> cell_cycle_time_initial_72_daughter;
        vector<double> cell_cycle_time_initial_73;
        vector<double> cell_cycle_time_initial_73_daughter;
        vector<double> cell_cycle_time_initial_74;
        vector<double> cell_cycle_time_initial_74_daughter;
        vector<double> cell_cycle_time_initial_75;
        vector<double> cell_cycle_time_initial_75_daughter;
        vector<double> cell_cycle_time_initial_76;
        vector<double> cell_cycle_time_initial_76_daughter;
        vector<double> cell_cycle_time_initial_77;
        vector<double> cell_cycle_time_initial_77_daughter;
        vector<double> cell_cycle_time_initial_78;
        vector<double> cell_cycle_time_initial_78_daughter;
        vector<double> cell_cycle_time_initial_79;
        vector<double> cell_cycle_time_initial_79_daughter;
        vector<double> cell_cycle_time_initial_80;
        vector<double> cell_cycle_time_initial_80_daughter;
        vector<double> cell_cycle_time_initial_81;
        vector<double> cell_cycle_time_initial_81_daughter;
        vector<double> cell_cycle_time_initial_82;
        vector<double> cell_cycle_time_initial_82_daughter;
        vector<double> cell_cycle_time_initial_83;
        vector<double> cell_cycle_time_initial_83_daughter;
        vector<double> cell_cycle_time_initial_84;
        vector<double> cell_cycle_time_initial_84_daughter;
        vector<double> cell_cycle_time_initial_85;
        vector<double> cell_cycle_time_initial_85_daughter;
        vector<double> cell_cycle_time_initial_86;
        vector<double> cell_cycle_time_initial_86_daughter;
        vector<double> cell_cycle_time_initial_87;
        vector<double> cell_cycle_time_initial_87_daughter;
        vector<double> cell_cycle_time_initial_88;
        vector<double> cell_cycle_time_initial_88_daughter;
        vector<double> cell_cycle_time_initial_89;
        vector<double> cell_cycle_time_initial_89_daughter;
        vector<double> cell_cycle_time_initial_90;
        vector<double> cell_cycle_time_initial_90_daughter;
        vector<double> cell_cycle_time_initial_91;
        vector<double> cell_cycle_time_initial_91_daughter;
        vector<double> cell_cycle_time_initial_92;
        vector<double> cell_cycle_time_initial_92_daughter;
        vector<double> cell_cycle_time_initial_93;
        vector<double> cell_cycle_time_initial_93_daughter;
        vector<double> cell_cycle_time_initial_94;
        vector<double> cell_cycle_time_initial_94_daughter;
        vector<double> cell_cycle_time_initial_95;
        vector<double> cell_cycle_time_initial_95_daughter;
        vector<double> cell_cycle_time_initial_96;
        vector<double> cell_cycle_time_initial_96_daughter;
        vector<double> cell_cycle_time_initial_97;
        vector<double> cell_cycle_time_initial_97_daughter;
        vector<double> cell_cycle_time_initial_98;
        vector<double> cell_cycle_time_initial_98_daughter;
        vector<double> cell_cycle_time_initial_99;
        vector<double> cell_cycle_time_initial_99_daughter;
        vector<double> cell_cycle_time_initial_100;
        vector<double> cell_cycle_time_initial_100_daughter;
        vector<double> cell_cycle_time_initial_101;
        vector<double> cell_cycle_time_initial_101_daughter;
        vector<double> cell_cycle_time_initial_102;
        vector<double> cell_cycle_time_initial_102_daughter;
        vector<double> cell_cycle_time_initial_103;
        vector<double> cell_cycle_time_initial_103_daughter;
        vector<double> cell_cycle_time_initial_104;
        vector<double> cell_cycle_time_initial_104_daughter;
        vector<double> cell_cycle_time_initial_105;
        vector<double> cell_cycle_time_initial_105_daughter;
        vector<double> cell_cycle_time_initial_106;
        vector<double> cell_cycle_time_initial_106_daughter;
        vector<double> cell_cycle_time_initial_107;
        vector<double> cell_cycle_time_initial_107_daughter;
        vector<double> cell_cycle_time_initial_108;
        vector<double> cell_cycle_time_initial_108_daughter;
        vector<double> cell_cycle_time_initial_109;
        vector<double> cell_cycle_time_initial_109_daughter;
        vector<double> cell_cycle_time_initial_110;
        vector<double> cell_cycle_time_initial_110_daughter;
        vector<double> cell_cycle_time_initial_111;
        vector<double> cell_cycle_time_initial_111_daughter;
        vector<double> cell_cycle_time_initial_112;
        vector<double> cell_cycle_time_initial_112_daughter;
        vector<double> cell_cycle_time_initial_113;
        vector<double> cell_cycle_time_initial_113_daughter;
        vector<double> cell_cycle_time_initial_114;
        vector<double> cell_cycle_time_initial_114_daughter;
        vector<double> cell_cycle_time_initial_115;
        vector<double> cell_cycle_time_initial_115_daughter;
        vector<double> cell_cycle_time_initial_116;
        vector<double> cell_cycle_time_initial_116_daughter;
        vector<double> cell_cycle_time_initial_117;
        vector<double> cell_cycle_time_initial_117_daughter;
        vector<double> cell_cycle_time_initial_118;
        vector<double> cell_cycle_time_initial_118_daughter;
        vector<double> cell_cycle_time_initial_119;
        vector<double> cell_cycle_time_initial_119_daughter;
        vector<double> cell_cycle_time_initial_120;
        vector<double> cell_cycle_time_initial_120_daughter;
        vector<double> cell_cycle_time_initial_121;
        vector<double> cell_cycle_time_initial_121_daughter;
        vector<double> cell_cycle_time_initial_122;
        vector<double> cell_cycle_time_initial_122_daughter;
        vector<double> cell_cycle_time_initial_123;
        vector<double> cell_cycle_time_initial_123_daughter;
        vector<double> cell_cycle_time_initial_124;
        vector<double> cell_cycle_time_initial_124_daughter;
        vector<double> cell_cycle_time_initial_125;
        vector<double> cell_cycle_time_initial_125_daughter;
        vector<double> cell_cycle_time_initial_126;
        vector<double> cell_cycle_time_initial_126_daughter;
        vector<double> cell_cycle_time_initial_127;
        vector<double> cell_cycle_time_initial_127_daughter;
        vector<double> cell_cycle_time_initial_128;
        vector<double> cell_cycle_time_initial_128_daughter;
        vector<double> cell_cycle_time_initial_129;
        vector<double> cell_cycle_time_initial_129_daughter;
        vector<double> cell_cycle_time_initial_130;
        vector<double> cell_cycle_time_initial_130_daughter;
        vector<double> cell_cycle_time_initial_131;
        vector<double> cell_cycle_time_initial_131_daughter;
        vector<double> cell_cycle_time_initial_132;
        vector<double> cell_cycle_time_initial_132_daughter;
        vector<double> cell_cycle_time_initial_133;
        vector<double> cell_cycle_time_initial_133_daughter;
        vector<double> cell_cycle_time_initial_134;
        vector<double> cell_cycle_time_initial_134_daughter;
        vector<double> cell_cycle_time_initial_135;
        vector<double> cell_cycle_time_initial_135_daughter;
        vector<double> cell_cycle_time_initial_136;
        vector<double> cell_cycle_time_initial_136_daughter;
        vector<double> cell_cycle_time_initial_137;
        vector<double> cell_cycle_time_initial_137_daughter;
        vector<double> cell_cycle_time_initial_138;
        vector<double> cell_cycle_time_initial_138_daughter;
        vector<double> cell_cycle_time_initial_139;
        vector<double> cell_cycle_time_initial_139_daughter;
        vector<double> cell_cycle_time_initial_140;
        vector<double> cell_cycle_time_initial_140_daughter;
        vector<double> cell_cycle_time_initial_141;
        vector<double> cell_cycle_time_initial_141_daughter;
        vector<double> cell_cycle_time_initial_142;
        vector<double> cell_cycle_time_initial_142_daughter;
        vector<double> cell_cycle_time_initial_143;
        vector<double> cell_cycle_time_initial_143_daughter;
        vector<double> cell_cycle_time_initial_144;
        vector<double> cell_cycle_time_initial_144_daughter;
        vector<double> cell_cycle_time_initial_145;
        vector<double> cell_cycle_time_initial_145_daughter;
        vector<double> cell_cycle_time_initial_146;
        vector<double> cell_cycle_time_initial_146_daughter;
        vector<double> cell_cycle_time_initial_147;
        vector<double> cell_cycle_time_initial_147_daughter;
        vector<double> cell_cycle_time_initial_148;
        vector<double> cell_cycle_time_initial_148_daughter;
        vector<double> cell_cycle_time_initial_149;
        vector<double> cell_cycle_time_initial_149_daughter;
        vector<double> cell_cycle_time_initial_150;
        vector<double> cell_cycle_time_initial_150_daughter;
        vector<double> cell_cycle_time_initial_151;
        vector<double> cell_cycle_time_initial_151_daughter;
        vector<double> cell_cycle_time_initial_152;
        vector<double> cell_cycle_time_initial_152_daughter;
        vector<double> cell_cycle_time_initial_153;
        vector<double> cell_cycle_time_initial_153_daughter;
        vector<double> cell_cycle_time_initial_154;
        vector<double> cell_cycle_time_initial_154_daughter;
        vector<double> cell_cycle_time_initial_155;
        vector<double> cell_cycle_time_initial_155_daughter;
        vector<double> cell_cycle_time_initial_156;
        vector<double> cell_cycle_time_initial_156_daughter;
        vector<double> cell_cycle_time_initial_157;
        vector<double> cell_cycle_time_initial_157_daughter;
        vector<double> cell_cycle_time_initial_158;
        vector<double> cell_cycle_time_initial_158_daughter;
        vector<double> cell_cycle_time_initial_159;
        vector<double> cell_cycle_time_initial_159_daughter;
        vector<double> cell_cycle_time_initial_160;
        vector<double> cell_cycle_time_initial_160_daughter;
        vector<double> cell_cycle_time_initial_161;
        vector<double> cell_cycle_time_initial_161_daughter;
        vector<double> cell_cycle_time_initial_162;
        vector<double> cell_cycle_time_initial_162_daughter;
        vector<double> cell_cycle_time_initial_163;
        vector<double> cell_cycle_time_initial_163_daughter;
        vector<double> cell_cycle_time_initial_164;
        vector<double> cell_cycle_time_initial_164_daughter;
        vector<double> cell_cycle_time_initial_165;
        vector<double> cell_cycle_time_initial_165_daughter;
        vector<double> cell_cycle_time_initial_166;
        vector<double> cell_cycle_time_initial_166_daughter;
        vector<double> cell_cycle_time_initial_167;
        vector<double> cell_cycle_time_initial_167_daughter;
        vector<double> cell_cycle_time_initial_168;
        vector<double> cell_cycle_time_initial_168_daughter;
        vector<double> cell_cycle_time_initial_169;
        vector<double> cell_cycle_time_initial_169_daughter;
        vector<double> cell_cycle_time_initial_170;
        vector<double> cell_cycle_time_initial_170_daughter;
    } data;

    // Random number distributions
    struct Distributions {
        std::normal_distribution<double> radius_dist;
        std::normal_distribution<double> nu_dist;
        std::normal_distribution<double> E_dist;
        std::normal_distribution<double> c_dist;
        std::normal_distribution<double> l_dist;
        std::normal_distribution<double> growth_rate_dist;
    } dists;

public:
    // Constructor
    CellGrowthSimulation() : gen(seed) {
        initializeDistributions();
    }

    // Initialize random number distributions
    void initializeDistributions() {
        dists.radius_dist = std::normal_distribution<double>(params.R_mean, params.R_stdv);
        dists.nu_dist = std::normal_distribution<double>(params.nu_mean, params.nu_stdv);
        dists.E_dist = std::normal_distribution<double>(params.E_mean, params.E_stdv);
        dists.c_dist = std::normal_distribution<double>(params.c_mean, params.c_stdv);
        dists.l_dist = std::normal_distribution<double>(params.l_mean, params.l_stdv);
        dists.growth_rate_dist = std::normal_distribution<double>(0, Constants::GR_RATE_RD_STDV);
    }

    // Core functionality methods
    void initializeSystem() {
        data_container_initialize_1();
        initializeParticles();
        initialize_position();
    }

    void runSimulation() {
        for (int t = params.startTime; t <= params.stopTime; t += params.timeStep) {
            calculate_force();
            time_advance();
            radius_change();
            
            if (t % params.timeStep == 0) {
                SnapshotDataPrint();
                PositionDataPrint();
            }
        }
    }

    // Periodic boundary conditions
    double PBC(double X) const {
        if (X > params.boxDimension / 2.0) {
            return X - params.boxDimension;
        } else if (X < -params.boxDimension / 2.0) {
            return X + params.boxDimension;
        }
        return X;
    }

    // Force calculation
    void calculate_force() {
        // Reset forces
        std::fill(data.force_x.begin(), data.force_x.end(), 0.0);
        std::fill(data.force_y.begin(), data.force_y.end(), 0.0);
        std::fill(data.force_z.begin(), data.force_z.end(), 0.0);

        // Calculate forces between particles
        for (int i = 0; i < state.numParticles; i++) {
            for (int j = i + 1; j < state.numParticles; j++) {
                calculate_pair_force(i, j);
            }
        }
    }

    // Calculate force between a pair of particles
    void calculate_pair_force(int i, int j) {
        state.dx = PBC(data.pos_x[i] - data.pos_x[j]);
        state.dy = PBC(data.pos_y[i] - data.pos_y[j]);
        state.dz = PBC(data.pos_z[i] - data.pos_z[j]);
        
        state.DIST_sq = state.dx * state.dx + state.dy * state.dy + state.dz * state.dz;
        state.DIST = std::sqrt(state.DIST_sq);
        
        double R_sum = data.radius_list[i] + data.radius_list[j];
        
        if (state.DIST < R_sum) {
            state.overlap = R_sum - state.DIST;
            
            // Calculate repulsive force
            double F_rep = calculate_repulsive_force(i, j);
            
            // Calculate adhesive force
            double F_ad = calculate_adhesive_force(i, j);
            
            // Total force
            double F_total = F_rep + F_ad;
            
            // Update forces
            double F_x = F_total * state.dx / state.DIST;
            double F_y = F_total * state.dy / state.DIST;
            double F_z = F_total * state.dz / state.DIST;
            
            data.force_x[i] += F_x;
            data.force_y[i] += F_y;
            data.force_z[i] += F_z;
            
            data.force_x[j] -= F_x;
            data.force_y[j] -= F_y;
            data.force_z[j] -= F_z;
        }
    }

    // Calculate repulsive force between particles
    double calculate_repulsive_force(int i, int j) const {
        double E_eff = 2.0 * data.elastic_moduli[i] * data.elastic_moduli[j] / 
                      (data.elastic_moduli[i] + data.elastic_moduli[j]);
        double nu_eff = 2.0 * data.poisson_ratio[i] * data.poisson_ratio[j] / 
                       (data.poisson_ratio[i] + data.poisson_ratio[j]);
        
        return (4.0/3.0) * E_eff * std::sqrt(data.radius_list[i] * data.radius_list[j] / 
               (data.radius_list[i] + data.radius_list[j])) * 
               std::pow(state.overlap, 1.5) * (1.0 - nu_eff);
    }

    // Calculate adhesive force between particles
    double calculate_adhesive_force(int i, int j) const {
        double c_factor = data.receptor_concentration[i] * data.ligand_concentration[j] + 
                         data.receptor_concentration[j] * data.ligand_concentration[i];
        return params.fad * c_factor * state.overlap;
    }

    // Time advancement
    void time_advance() {
        for (int i = 0; i < state.numParticles; i++) {
            // Update velocities
            data.vel_x[i] += params.kdt * data.force_x[i];
            data.vel_y[i] += params.kdt * data.force_y[i];
            data.vel_z[i] += params.kdt * data.force_z[i];
            
            // Update positions
            data.pos_x[i] += params.kdt * data.vel_x[i];
            data.pos_y[i] += params.kdt * data.vel_y[i];
            data.pos_z[i] += params.kdt * data.vel_z[i];
            
            // Apply periodic boundary conditions
            data.pos_x[i] = PBC(data.pos_x[i]);
            data.pos_y[i] = PBC(data.pos_y[i]);
            data.pos_z[i] = PBC(data.pos_z[i]);
        }
    }

    // Initialization methods
    void data_container_initialize_1() {
        data.radius_list.resize(params.N0);
        data.pos_x.resize(params.N0);
        data.pos_y.resize(params.N0);
        data.pos_z.resize(params.N0);
        data.vel_x.resize(params.N0);
        data.vel_y.resize(params.N0);
        data.vel_z.resize(params.N0);
        data.force_x.resize(params.N0);
        data.force_y.resize(params.N0);
        data.force_z.resize(params.N0);
        data.elastic_moduli.resize(params.N0);
        data.poisson_ratio.resize(params.N0);
        data.receptor_concentration.resize(params.N0);
        data.ligand_concentration.resize(params.N0);
        data.growth_rate.resize(params.N0);
        data.growth_rate_daughter.resize(params.N0);
        data.cell_cycle_time.resize(params.N0);
        data.cell_cycle_time_daughter.resize(params.N0);
        data.cell_cycle_time_initial.resize(params.N0);
        data.cell_cycle_time_initial_daughter.resize(params.N0);
    }

    void initializeParticles() {
        for (int i = 0; i < params.N0; i++) {
            // Initialize particle properties
            data.radius_list[i] = dists.radius_dist(gen);
            data.elastic_moduli[i] = dists.E_dist(gen);
            data.poisson_ratio[i] = dists.nu_dist(gen);
            data.receptor_concentration[i] = dists.c_dist(gen);
            data.ligand_concentration[i] = dists.c_dist(gen);
            
            // Initialize growth rates
            data.growth_rate[i] = params.rv + dists.growth_rate_dist(gen);
            data.growth_rate_daughter[i] = params.rv + dists.growth_rate_dist(gen);
            
            // Initialize cell cycle times
            data.cell_cycle_time[i] = 0.0;
            data.cell_cycle_time_daughter[i] = 0.0;
            data.cell_cycle_time_initial[i] = 0.0;
            data.cell_cycle_time_initial_daughter[i] = 0.0;
        }
    }

    void initialize_position() {
        // Initialize positions randomly within the box
        std::uniform_real_distribution<double> pos_dist(-params.boxDimension/2.0, params.boxDimension/2.0);
        
        for (int i = 0; i < params.N0; i++) {
            data.pos_x[i] = pos_dist(gen);
            data.pos_y[i] = pos_dist(gen);
            data.pos_z[i] = pos_dist(gen);
            
            // Initialize velocities to zero
            data.vel_x[i] = 0.0;
            data.vel_y[i] = 0.0;
            data.vel_z[i] = 0.0;
        }
    }

    // Data output methods
    void SnapshotDataPrint() const {
        string filename = "snapshot_data_" + to_string(params.startTime) + ".dat";
        ofstream outfile(filename);
        
        if (!outfile) {
            throw runtime_error("Failed to open snapshot data file for writing");
        }
        
        outfile << fixed << setprecision(6);
        for (int i = 0; i < state.numParticles; i++) {
            outfile << setw(10) << data.pos_x[i] << " "
                   << setw(10) << data.pos_y[i] << " "
                   << setw(10) << data.pos_z[i] << " "
                   << setw(10) << data.radius_list[i] << endl;
        }
    }

    void PositionDataPrint() const {
        string filename = "position_data_" + to_string(params.startTime) + ".dat";
        ofstream outfile(filename);
        
        if (!outfile) {
            throw runtime_error("Failed to open position data file for writing");
        }
        
        outfile << fixed << setprecision(6);
        for (int i = 0; i < state.numParticles; i++) {
            outfile << setw(10) << data.pos_x[i] << " "
                   << setw(10) << data.pos_y[i] << " "
                   << setw(10) << data.pos_z[i] << endl;
        }
    }

    // Cleanup methods
    void cleanup() {
        data.radius_list.clear();
        data.pos_x.clear();
        data.pos_y.clear();
        data.pos_z.clear();
        data.vel_x.clear();
        data.vel_y.clear();
        data.vel_z.clear();
        data.force_x.clear();
        data.force_y.clear();
        data.force_z.clear();
        data.elastic_moduli.clear();
        data.poisson_ratio.clear();
        data.receptor_concentration.clear();
        data.ligand_concentration.clear();
        data.growth_rate.clear();
        data.growth_rate_daughter.clear();
        data.cell_cycle_time.clear();
        data.cell_cycle_time_daughter.clear();
        data.cell_cycle_time_initial.clear();
        data.cell_cycle_time_initial_daughter.clear();
    }

    // Main function
    static int main(int argc, char *argv[]) {
        try {
            CellGrowthSimulation simulation;
            simulation.initializeSystem();
            simulation.runSimulation();
            simulation.cleanup();
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
};

// Global main function
int main(int argc, char *argv[]) {
    return CellGrowthSimulation::main(argc, argv);
}
