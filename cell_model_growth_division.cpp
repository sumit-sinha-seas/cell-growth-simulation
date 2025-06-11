#include <bits/stdc++.h>
#include <chrono>
// #include"connection_number.h"

using namespace std;

using namespace std::chrono;

random_device r;
seed_seq seed{r(), r(), r(), r(), r()};
// This defines the seed that we are using for generating random number.
mt19937 gen(seed);

/*CONSTANT PARAMETERS*/
/*************************************************************************************************************************************/
int N0 = 0;
// N0 denotes the initial number of particles.

double boxDimension = 0.0;
// boxDimension defines the dimension of the cubic box.

double fad = 0;
// fad denotes the cell-cell adhesion strength.

int sp_dim = 3;

int startTime = 0;

int stopTime = 0;

int time_step = 10;

int SC_startTime = 10;

int DynamicVectorDim = 0;

int id_lv = 1, id_DV = 2, id_dd = 3;

// int minDataPoints = 140;

int minDataPoints = 170;

double eta = 5e-3;
// eta denoes ECM viscosity.

double l_mean=0;
double l_stdv=25;

double R_mean = 4.5;
double R_stdv = 0.5;
// R_mean & R_stdv are the parameters that are used to build the normal distribution which will generate the radii list initially.

double mu_mean = 0.045;
// mu is self propulsion mobility.

double nu_mean = 0.5;
double nu_stdv = 0.02;
// nu is poisson ratio. nu_mean & nu_stdv are the parameters that are used to calculate the repulsive force. These represent properties of a particular cell.
// These are updated at each cell cycle.

double E_mean = 1e-3;
double E_stdv = 1e-4;
// E is elastic moduli. E_mean & E_stdv are the parameters that are used to calculate the repulsive force. These represent propeties of a particular cell.
// These are updated at each cell cycle.

double c_mean = 0.9;

// double c_mean=1.0;

double c_stdv = 0.02;
// c1 denotes the normalized receptor concentration, c2 denotes the normalized ligand concentration. c1 & c2 values are taken from a uniform distribution that is
//  defined by c_mean & c_stdv. These represents properties of a particular cell. c1 & c2 are updated at each cell cycle.

double gamma_max = 1e-4;
// gamma_max denotes ahesive friction.

double Rm = 5.0;
// Rm denotes mitotic radius.

double CC_time = 54000;
//

double Pc = 1e-4;
// double Pc = 1e10;
//   Pc denotes critical pressure.

double grrateRd_stdv = 1e-5;
// grrated_stdv denotes the standard deviation of the growth rate distribution.

// double rv = (2 * M_PI * pow(Rm, 3)) / (3 * CC_time);

double rv = (2 * M_PI * pow(Rm, 3)) / (3 * CC_time);
// rv is a constant which appears in the expression of mean growth rate.

double eps = 1e-30;

double sigma = 0.01;

double kdt = ((double)time_step) * (1e-6);

double scaleFactor = 0.0;

double flag_1 = -50.0;

double flag_2 = -10.0;

// double R_dot=1/(4*M_PI);

/*****************************************************************************/

/*VARIABLE PARAMETERS*/
/**************************************************************************************************************/
double phi = 0; // phi defines the packing fraction.
// double sigma = 0.7;   // sigma defines the extent upto which two cells can overlap initially
int numParticles = 0; //  numParticles defines the total number of cells at each time step.
// double boxDimension = 50.0; // boxDimension defines the dimension of the cubic box.

int dd_numParticles = 0;

int DV_numParticles = 0;

double dx = 0, dy = 0, dz = 0, DIST = 0, DIST_sq = 0;
// dx can denote two things. (a) dx=x_i-x_j (b) dx=x_i(t2)-x_i(t1)
//  DISt can denote two things. (a) euclidean distance between two particles. (b) euclidean distance between the position vector of particle at time t1 &
//  the postion vetor of the same particle at time t2
// DIST_sq denotes the sqrt of DIST
double R_cutoff = 0.0;
// This defines the cutoff distance to make the code more efficient.

double Rd = 0.0;
// This defines the radius of each particle.

double delta_Rd = 0.0;
// This denotes the change in radius of each particle.

double Rdtr = 0.0;
// Rdtr denotes the radius of the daughter particle.

double PRESSURE = 0.0;
// It stores the pressure value temporarily.

double E = 0.0, nu = 0.0, c1 = 0.0, c2 = 0.0;
// These variable store elastic moduli, poisson ratio, receptor & ligand concentration temporarily.

double max_DIA = 0.0;
// max_DIA denotes the maximum diameter at current time.

// double rds1=0.0, rds2=0.0;

// double inv_num_density=0;

/*****************************************************************************************************************/

/*Correlation Function Parameters*/

// double corr_factor_1 = 0;
/*double corr_up_limit = 400.0;
double corr_low_limit = 0.0;
double corr_increment = 0.01;
int corr_step = (int)((corr_up_limit - corr_low_limit) / corr_increment);*/

/*********************************************************************/

/*double pd = 0.0, pd1 = 0.0, pd2 = 0.0;
double par1_x = 0, par2_x = 0, par1_y = 0, par2_y = 0;*/

double overlap = 0.0;

// double overlap = 0.0;
//  overlap represents the overlap length between two cells.

// double cfactor_total=0;
// It represents the total c_factor for each particle.

// pd1 defines the (x1)-(x2) for two particles 1 & 2.
// pd2 defines the (y1)-(y2) for two particles 1 & 2.
// double pair_distance = 0.0;
// pd defines the euclidean distance between two particles 1 & 2.

/**CLA VARIABLES*/

int ensemblesize = 1;
int filenum = 1;

/*************************************************/

// double h_ij = 0;
// double R1 = 0, R2 = 0;
// double F_factor2 = 0, F_factor3 = 0;

// vector<double> Distance;
vector<double> radius_list;
// radius_list data container stores all the radii value.

vector<double> pos_x;
// pos_x data container stores the x-coordinates of the particles.

vector<double> pos_y;
// pos_y data container stores the y-coordinates of the particles.

vector<double> pos_z;
// pos_z data container stores the z-coordinates of the particles.

// vector<double> vel_x;
//  vel_x data container stores the x-component of the velocities of the particles at time t. It is initialized at every time step.

// vector<double> vel_y;
//  vel_y data container stores the y-component of the velocities of the particles at time t. It is initialized at every time step.

vector<double> force_x;
// force_x data container stores the x-component of the force on the particle i at time t. It is initialized at every time step.

vector<double> force_y;
// force_y data container stores the y-component of the force on the particle i at time t. It is initialized at every time step.

vector<double> force_z;
// force_z data container stores the z-component of the force on the particle i at time t. It is initialized at every time step.

vector<double> pressure;
// Pressure vector stores the pressure value for particle i. It is initialized at every time step.

/*
############

vector<double> pos_x0;
// pos_x0 data container stores the x0-coordinates of the particles.

vector<double> pos_y0;
// pos_y0 data container stores the y0-coordinates of the particles.

vector<double> pos_z0;
// pos_z0 data container stores the z0-coordinates of the particles.

*/

vector<double> gammaMt_vsc;
// gammaMt_vsc stores the values that represent friction due to viscosity of the medium.

vector<double> gammaMt_ad;
// gammaMt_ad

vector<double> gammaMt_adx1;
// the i-th element of gammaMt_adx1 represents the sum of Aij*c1(i)*c2(j)+c2(i)*c1(j)*the component of the unit vector in the x-direction.

vector<double> gammaMt_ady1;
// the i-th element of gammaMt_adx1 represents the sum of Aij*c1(i)*c2(j)+c2(i)*c1(j)*the component of the unit vector in the y-direction.

vector<double> gammaMt_adz1;
// the i-th element of gammaMt_adx1 represents the sum of Aij*c1(i)*c2(j)+c2(i)*c1(j)*the component of the unit vector in the z-direction.

vector<double> gammaMt_ad2;

vector<double> nu_list;

vector<double> E_list;

vector<double> c1_list;

vector<double> c2_list;

vector<double> uv;

vector<int> particle_id;

/*set<int> C1;

set<int> C2;

set<int> C3;*/

// vector<vector<int>> nn_list;
//  nn_list is a 2D data container that stores the neighbor a particular cell at each time step.

vector<double> TIME;

// vector<int> timePRINT_indx;

// set<int> timePRINT_indxSet;

vector<int> timePRINT_chck;

//vector<double> PARTICLE_NUM;
vector<int> PARTICLE_NUM;

/*

###########

//vector<double> MSD_TIME;

//vector<double> CORR_R;

//vector<double> CORR_FUNC;

*/

normal_distribution<double> dist1(R_mean, R_stdv);

normal_distribution<double> dist2(nu_mean, nu_stdv);

normal_distribution<double> dist3(E_mean, E_stdv);

normal_distribution<double> dist4(c_mean, c_stdv);

normal_distribution<double> dist5(l_mean, l_stdv);

uniform_real_distribution<double> dist6(-1, 1);

uniform_real_distribution<double> dist12(0, 1);

// gamma_not is friction coefficient

// double F_factor1 = (2 * E) / (3 * (1 - (nu * nu)));
/*double tau = (mean_R * mean_R) / (mu * mu);

int itr_steps = (int)tau;*/

// double particleRadius=0.0;
// void data_container_initialize_1()
//{}

void data_container_initialize_1()
{

    particle_id.assign(N0, id_lv);

    force_x.assign(N0, 0.0);

    force_y.assign(N0, 0.0);

    force_z.assign(N0, 0.0);

    pos_x.assign(N0, 0.0);

    pos_y.assign(N0, 0.0);

    pos_z.assign(N0, 0.0);

    pressure.assign(N0, 0.0);

    radius_list.assign(N0, 0.0);

    gammaMt_vsc.assign(N0, 0.0);

    nu_list.assign(N0, 0.0);

    E_list.assign(N0, 0.0);

    c1_list.assign(N0, 0.0);

    c2_list.assign(N0, 0.0);

    gammaMt_adx1.assign(N0, 0.0);

    gammaMt_ady1.assign(N0, 0.0);

    gammaMt_adz1.assign(N0, 0.0);

    gammaMt_ad2.assign(N0, 0.0);

    // vel_x.resize(numParticles, 0.0);

    // vel_y.resize(numParticles, 0.0);

    // cout << "force_x_size:" << force_x.size() << endl;
}

void DynamicDataSavings1()
{
    // int DynamicVectorDim = ((stopTime - startTime) / time_step) + 1;
    int CCT=(int)(CC_time/10.0);

    cout<<CCT<<endl;

    for (int i = 0; i < DynamicVectorDim; i++)
    {
        TIME.push_back((double)(i * time_step));
        timePRINT_chck.push_back(-1);
        // PARTICLE_NUM.push_back(0);
        //  MSD_TIME.push_back(0);
    }

    for (int i = 0; i < DynamicVectorDim; i++)
    {
        if(i%CCT==0)
            timePRINT_chck[i]=1;
        // PARTICLE_NUM.push_back(0);
        //  MSD_TIME.push_back(0);
    }

    int SC_stopTime = stopTime / time_step;
    int value = SC_startTime;
    double scaleFactor = pow((double)(SC_stopTime / SC_startTime), 1.0 / (minDataPoints - 1));
    // cout<<value<<" "<<SC_stopTime<<" "<<scaleFactor<<endl;
    while (value <= SC_stopTime)
    {
        // timePRINT_indxSet.insert(value);
        // cout<<"dd container works"<<endl;
        timePRINT_chck[value] = 1;
        value = (round(value * scaleFactor));
        // timePRINT_chck[value] = 1;
    }

    // for (const int &i : timePRINT_indxSet)
    // timePRINT_indx.push_back(i);
}

void DynamicDataSavings2()
{
    for (int i = 0; i < DynamicVectorDim; i++)
    {

        PARTICLE_NUM.push_back(0);
        // MSD_TIME.push_back(0);
    }
}

/*void data_container_initialize_2()
{
    for (int k = 0; k < itr_steps; k++)

    {

        MSD_time.push_back(0.0);

        TIME.push_back((double)(k * time_step));
    }

    cout << "MSD_time_size:" << MSD_time.size() << endl;
}

void data_container_initialize_3()
{
    for (int k = 0; k < corr_step; k++)

    {

        CORR_func.push_back(0.0);

        CORR_R.push_back((double)k * (double)corr_increment);
    }

    cout << "CORR_func_size:" << CORR_func.size() << endl;
}*/

double PBC(double X)
{

    // double x=0;
    if (X > boxDimension)

    {
        while (X > boxDimension)

            X = X - boxDimension;

        return X;
    }

    else if (X < 0)
    {
        while (X < 0)

            X = X + boxDimension;

        return X;
    }

    else

        return X;
}

/*

double min_dist(double x1, double x2)
{

    double dx = x1 - x2;

    dx = dx - round(dx / boxDimension) * dx;

    return dx;


}*/

/*double pair_distance(double x1, double x2, double y1, double y2)
{

    // p1=

    double d1 = min_dist(x1, x2);

    double d2 = min_dist(y1, y2);

    double D = sqrt((d1 * d1) + (d2 * d2));

    return D;
}*/

/*double normal1(double mean, double stdv)
{

    double temp;

    normal_distribution<double> dist1(mean, stdv);

    while (true)
    {
        temp = dist1(gen);

        if (temp > 0)

            break;
    }

    return temp;
}

double normal2(double mean, double stdv)
{

    double temp;

    normal_distribution<double> dist2(mean, stdv);

    temp = dist2(gen);

    return temp;
}

double normal3(double mean, double stdv)
{

    double temp;

    normal_distribution<double> dist3(mean, stdv);

    while (true)
    {
        temp = dist3(gen);

        if (temp >= eps && temp <= 1)

            break;
    }

    temp = dist3(gen);

    return temp;
}*/

/*double randouble(double min, double max)
{

    double temp;

    uniform_real_distribution<double> dist4(min, max);

    temp = dist4(gen);

    return temp;
}*/

void unit_vector()
{

    uv.assign(sp_dim, 0.0);

    double x, y, z;
    do
    {
        x = dist6(gen);
        y = dist6(gen);
        z = dist6(gen);
        DIST_sq = x * x + y * y + z * z;
    } while (DIST_sq > 1.0);

    uv[0] = x / sqrt(DIST_sq);

    uv[1] = y / sqrt(DIST_sq);

    uv[2] = z / sqrt(DIST_sq);
}

void initializeParticles()
{

    double temp_R, temp_nu, temp_E, temp_c1, temp_c2;

    //double sum_radi_sqr = 0;

    /*normal_distribution<double> dist1(R_mean, R_stdv);

    normal_distribution<double> dist2(nu_mean, nu_stdv);

    normal_distribution<double> dist3(E_mean, E_stdv);

    normal_distribution<double> dist4(c_mean, c_stdv);*/

    for (int i = 0; i < N0; i++)

    {

        while (true)
        {
            temp_R = dist1(gen);

            // cout<<temp_R<<endl;

            // if (temp_R > 0 && temp_R<=5.0)
            if (temp_R > 0)

                break;
        }

        radius_list[i] = temp_R;

        //sum_radi_sqr += temp_R * temp_R;

        gammaMt_vsc[i] = eta * temp_R * 6 * M_PI;
    }

    // cout<<"R, eta okay"<<endl;

    for (int i = 0; i < N0; i++)

    {
        nu_list[i] = dist2(gen);
    }

    // cout<<"nu okay"<<endl;

    for (int i = 0; i < N0; i++)

    {

        while (true)

        {

            temp_E = dist3(gen);

            if (temp_E > 0)

                break;
        }

        E_list[i] = temp_E;
    }

    // cout<<"E ok"<<endl;

    for (int i = 0; i < N0; i++)
    {
        while (true)

        {

            temp_c1 = dist4(gen);

            if (temp_c1 >= 0 && temp_c1 <= 1)

                break;
        }

        c1_list[i] = temp_c1;
    }

    for (int i = 0; i < N0; i++)
    {
        while (true)

        {

            temp_c2 = dist4(gen);

            if (temp_c2 >= 0 && temp_c2 <= 1)

                break;
        }

        c2_list[i] = temp_c2;
    }

    // cout<<"c1,c2 ok"<<endl;

    //phi = (M_PI * sum_radi_sqr) / (pow(boxDimension, 3));
}

void initialize_position()
{

   for (int i = 0; i < N0; i++)
   {

        pos_x[i] = dist5(gen);

        pos_y[i] = dist5(gen);

        pos_z[i] = dist5(gen);

    }

    // cout << "inital position ok" << endl;
}

void calculate_force()
{

    double Aij, sqronebyR, onebyR, hij, cfactor_1;

    double force_rpl, force_ads, force_ij;

    double temp1, temp2, temp3, temp4, temp6;

    int temp_id = 0;

    double n1, n2, n3;

    force_x.assign(numParticles, 0.0);

    force_y.assign(numParticles, 0.0);

    force_z.assign(numParticles, 0.0);

    pressure.assign(numParticles, 0.0);

    gammaMt_adx1.assign(numParticles, 0.0);

    gammaMt_ady1.assign(numParticles, 0.0);

    gammaMt_adz1.assign(numParticles, 0.0);

    gammaMt_ad2.assign(numParticles, 0.0);

    /*if (force_x.size() != numParticles)
        cout << "error" << endl;*/

    /*for (int i = 0; i < numParticles; i++)
    {
        force_x[i]=0.0;
        force_y[i] = 0.0;
        force_z[i] = 0.0;
        pressure[i] = 0.0;
        gammaMt_adx1[i] = 0.0;
        gammaMt_ady1[i] = 0.0;
        gammaMt_adz1[i] = 0.0;
        gammaMt_ad2[i] = 0.0;
    }*/

    double rds1, rds2;

    for (int i = 0; i < numParticles - 1; i++)
    {

        double gamma_ad1 = 0, gamma_ad2 = 0;

        // cout<<i<<",";

        rds1 = radius_list[i];

        temp_id = particle_id[i];

        if (temp_id == id_dd)
        // if (rds1 < flag_2)
        {
            // cout << "Force Function should not work at i:" << i << endl;
            // C1.insert(i);
            continue;
        }

        else
        {
            for (int j = i + 1; j < numParticles; j++)

            {

                // cout<<j<<","<<endl;

                rds2 = radius_list[j];

                temp_id = particle_id[j];

                if (temp_id == id_dd)
                // if (rds2 < flag_2)
                {

                    // C1.insert(j);
                    //  cout << "Force Function should not work at j:" << j << endl;
                    continue;
                }

                else
                {
                    dx = pos_x[i] - pos_x[j];

                    dy = pos_y[i] - pos_y[j];

                    dz = pos_z[i] - pos_z[j];

                    /*dx = dx - round(dx / boxDimension) * boxDimension;
                    dy = dy - round(dy / boxDimension) * boxDimension;
                    dz = dx - round(dz / boxDimension) * boxDimension;*/

                    DIST = sqrt((dx * dx + dy * dy + dz * dz));

                    // sqronebyR = (1.0 / radius_list[i]) + (1.0 / radius_list[j]);

                    sqronebyR = (1.0 / rds1) + (1.0 / rds2);

                    // temp6 = radius_list[i] + radius_list[j];

                    temp6 = rds1 + rds2;

                    overlap = temp6 - DIST;

                    // if(overlap<0)
                    // cout<<"overlap:"<<overlap<<endl;

                    // cout << "overlap:" << overlap << endl;

                    if (overlap > 0.0)
                    {

                        onebyR = sqrt(sqronebyR);
                        hij = overlap;
                        temp1 = pow(hij, 1.5);
                        Aij = M_PI * hij * (1.0 / sqronebyR);
                        // sqrt_Aij = sqrt(Aijt);
                        temp2 = 0.75 * (((1.0 - nu_list[i] * nu_list[i]) / E_list[i]) + ((1.0 - nu_list[j] * nu_list[j]) / E_list[j])) * onebyR;
                        force_rpl = temp1 / temp2;
                        cfactor_1 = c1_list[i] * c2_list[j] + c1_list[j] * c2_list[i];
                        temp3 = Aij * cfactor_1;
                        force_ads = 0.5 * temp3 * fad;
                        force_ij = force_rpl - force_ads;
                        n1 = dx / DIST;
                        n2 = dy / DIST;
                        n3 = dz / DIST;
                        force_x[j] -= (force_ij)*n1;
                        force_x[i] += (force_ij)*n1;
                        force_y[j] -= (force_ij)*n2;
                        force_y[i] += (force_ij)*n2;
                        force_z[j] -= (force_ij)*n3;
                        force_z[i] += (force_ij)*n3;
                        pressure[i] += (abs(force_ij) / Aij);
                        pressure[j] += (abs(force_ij) / Aij);

                        gammaMt_adx1[i] += n1 * temp3;

                        gammaMt_adx1[j] -= n1 * temp3;

                        gammaMt_ady1[i] += n2 * temp3;

                        gammaMt_ady1[j] -= n2 * temp3;

                        gammaMt_adz1[i] += n3 * temp3;

                        gammaMt_adz1[j] -= n3 * temp3;

                        gammaMt_ad2[i] += temp3;
                        gammaMt_ad2[j] += temp3;
                    }
                }
            }

            /*if (j == numParticles - 1)
            {

                temp4 = force_x[i] * force_x[i] + force_y[i] * force_y[i] + force_z[i] * force_z[i];

                if (temp4 != 0)

                {
                    gamma_ad1 = gamma_ad1 / sqrt(temp4);

                    gamma_ad1 = gamma_ad1 + 1;

                    gammaMt_ad[i] = 0.25 * gamma_max * (gamma_ad1 + gamma_ad2);
                }

                else
                    gammaMt_ad[i] = 0;*/
        }
    }
}

/*double calc_force_factor_1(double H, double R1, double R2)
{

    double F_factor2 = sqrt((R1 * R2) / (R1 + R2));

    double F_factor3 = pow(H, 1.5) * F_factor2 * F_factor1;

    return F_factor3;
}*/

void initialize_system()
{

    numParticles = 0;

    dd_numParticles = 0;

    DV_numParticles = 0;

    data_container_initialize_1();

    initializeParticles();

    // cout << "initial 1 ok" << endl;

    initialize_position();

    // cout << "initial 2 ok" << endl;

    numParticles = N0;

    calculate_force();
}

// check the enext two functions later.

/*void inititalDataPrint()
{

    for (int i = 0; i < radius_list.size(); i++)

    {
        cout << "Radius:" << radius_list[i] << endl;

        // cout<<"overlap:"<<overlap<<endl;

        cout << "Fx:" << force_x[i] << endl;

        cout << "Fy:" << force_y[i] << endl;

        cout << "Fz:" << force_z[i] << endl;

        cout << "P:" << pressure[i] << endl;

        cout<<"id:"<< particle_id[i]<<endl;
    }
}*/

void SnapshotDataPrint()
{

    for (int i = 0; i < numParticles; i++)

    {
        cout << "Radius:" << radius_list[i] << endl;

        cout << "Fx:" << force_x[i] << endl;

        cout << "Fy:" << force_y[i] << endl;

        cout << "Fz:" << force_z[i] << endl;

        cout << "P:" << pressure[i] << endl;

        // cout << "P:" << pressure[i] << endl;

        cout << "id:" << particle_id[i] << endl;
    }
}

void PositionDataPrint()
{

    for (int i = 0; i < numParticles; i++)

    {
        cout << "rx:" << pos_x[i] << " "
        << "ry:" << pos_y[i] << " "
        << "rz:" << pos_z[i] << " " << endl;
    }
}

/*double calc_msd(double x1, double x2, double y1, double y2)
{

    // 1 denotes first particle of the pair, 2 denotes second particle of the pair.

    // x1 denotes the x coordinate of the first particle.
    // x2 denotes the x coordinate of the second particle.
    // y1 denotes the y coordinate of the first particle.
    // y2 denotes the y coordinate of the second particle.

    double MSD = pair_distance(x1, x2, y1, y2) * pair_distance(x1, x2, y1, y2);

    return MSD;
}*/

/*void correlation_factor1(double p, double n)
{

    for(int i=0;i<CORR_func.size();i++)
    {

     //double p1 = 0.0;

    if ((CORR_R[i] - p) > -1e-3 && (CORR_R[i] - p) < 1e-3)

        CORR_func[i]+=1*n;

    else

        continue;

    }

}*/

/*void calc_dynamics()

{

    double q1 = 0, q2 = 0, s1 = 0, s2 = 0, Rd1 = 0, Rd2 = 0, overlap = 0;

    double noise = 0.0;

    double inv_num_density = 0;

    int check_1 = 0;

    initializeParticles();

    for (int m = 0; m < numParticles; m++)

    {

        pos_x0.push_back(pos_x[m]);

        pos_y0.push_back(pos_y[m]);
    }

    cout << "pos_x0_size:" << pos_x0.size() << endl;

    inv_num_density = (boxDimension * boxDimension) / (double)(numParticles * numParticles);

    cout << "inv_num_density:" << inv_num_density << endl;

    // 1 denotes particle 1, 2 denotes particle 2.
    // q is equivalent to x
    // s is equivalent to y

    // q1 denotes the x coordinate of the first particle.
    // q2 denotes the x coordinate of the second particle.
    // s1 denotes the y coordinate of the first particle.
    // s2 denotes the y coordinate of the second particle.

    normal_distribution<double> dist3(0, 1);

    // int_itr_steps = 1;

    // double time_total=10*tau;

    for (int k = 0; k < itr_steps; k++)

    {

        // double noise = dist3(gen);

        // MSD_time[k].push_back(0);

        if (k == 0)

        {

            /*initializeParticles();

            for (int m = 0; m < numParticles; m++)

            {

                pos_x0.push_back(pos_x[m]);

                pos_y0.push_back(pos_y[m]);
            }

            cout << "pos_x0_size:" << pos_x0.size() << endl;

            inv_num_density = (boxDimension * boxDimension) / (double)(numParticles * numParticles);

            cout << "inv_num_density:" << inv_num_density << endl;

            continue;
        }
        // else if (k != 0 && k % 10 == 0)

        else

        {

            // vector<vector<int>>().swap(nn_list);

            for (int i = 0; i < numParticles - 1; ++i)

            {

                // Here, i denotes the first particle of the pair, & j denotes the second particle of the pair.

                // nn_list[i].clear();

                q1 = pos_x[i];

                // q1 denotes the x coordinate of the first particle.

                s1 = pos_y[i];

                // s1 denotes the y coordinate of the first particle.

                Rd1 = radius_list[i];

                // Rd1 denotes the radius of the first particle.

                for (int j = i + 1; j < numParticles; ++j)

                {

                    // vector<int>().swap(nn_list[j]);

                    q2 = pos_x[j];

                    // q2 denotes the x coordinate of the second particle.

                    s2 = pos_y[j];

                    // s2 denotes the y coordinate of the second particle.

                    Rd2 = radius_list[j];

                    // Rd2 denotes the radius of the second particle.

                    //
                    // correlation_factor1(pair_distance(q1, q2, s1, s2),inv_num_density);

                    pd = pair_distance(q1, q2, s1, s2);

                    overlap = (Rd1 + Rd2) - pd;

                    if (k == (int)(itr_steps - 10))
                    {

                        check_1++;

                        int binIndex = (int)(pd / corr_increment);

                        if (binIndex < corr_step)

                            CORR_func[binIndex] += inv_num_density;
                    }

                    // overlap = calc_force_factor_1(Rd1, Rd2, q1, q2, s1, s2);

                    if (overlap > 0)

                    {

                        // nn_list[i].push_back(j);

                        // nn_list[j].push_back(i);

                        double ff_2 = calc_force_factor_1(overlap, Rd1, Rd2);

                        double z1 = min_dist(q1, q2);

                        double z2 = min_dist(s1, s2);

                        double z3 = sqrt((z1 * z1) + (z2 * z2));

                        double z4 = ff_2 * (z1 / z3);

                        double z5 = ff_2 * (z2 / z3);

                        force_x[i] += z4;

                        force_x[j] += -z4;

                        force_y[i] += z5;

                        force_y[j] += -z5;
                    }

                    else

                        continue;
                }

                noise = dist3(gen);

                vel_x[i] = force_x[i] / (radius_list[i] * gamma_not) + mu * noise;

                force_x[i] = 0;

                noise = dist3(gen);

                vel_y[i] = force_y[i] / (radius_list[i] * gamma_not) + mu * noise;

                force_y[i] = 0;

                double delta_x = vel_x[i] * time_step;

                double delta_y = vel_y[i] * time_step;

                /*if (k==1 && i<=5)

                {

                    cout << i<<" "<<"vx_i" <<vel_x[i]<< " " << "vy_i" << " " << vel_y[i] << endl;

                     //cout << "vx_2" <<vel_x[2]<< " " << "vy_2" << " " << vel_y[2] << endl;

                }*/

// double temp_x1=pos_x[i] + vel_x[i] * time_step;

/*if (delta_x > 0.5 * boxDimension)

    delta_x = delta_x - boxDimension;

else if (delta_x > -0.5 * boxDimension)

    delta_x = delta_x + boxDimension;

else

    delta_x = delta_x + 0.0;*/

/*pos_x[i] = PBC(pos_x[i] + delta_x);



pos_y[i] = PBC(pos_y[i] + delta_y);

MSD_time[k] += ((pos_x[i] - pos_x0[i]) * (pos_x[i] - pos_x0[i]) + (pos_y[i] - pos_y0[i]) * (pos_y[i] - pos_y0[i])) / (double)numParticles;*/

// MSD_time[k] += (pair_distance(pos_x[i], pos_x0[i], pos_y[i], pos_y0[i]) * pair_distance(pos_x[i], pos_x0[i], pos_y[i], pos_y0[i])) / (double)numParticles;
/*

            }
        }
    }

    cout << "check_1:" << check_1 << endl;
}*/

void radius_change()
{

    double grrateRd_mean;

    double temp, temp_DIST, temp_DIA, temp_Rd, temp2, temp3, temp4, temp5;
    double temp_check, temp_id;

    double temp_x, temp_y, temp_z;

    // double N0_t;

    int N0_t = numParticles;

    // cout << N0_t << endl;

    for (int i = 0; i < numParticles; i++)
    {

        Rd = radius_list[i];
        PRESSURE = pressure[i];
        temp_id = particle_id[i];

        if (temp_id == id_dd)
        // if (Rd < flag_2)
        {

            // cout << "id of dd particle:"<< temp_id << endl;

            // C2.insert(i);

            if (i == (N0_t - 1))

            {

                break;

                // continue;
            }

            else
                continue;
        }

        else
        {

            // if (Rd >= Rm && PRESSURE < Pc)

            // cout << i <<" "<< "Radius:" << Rd << "Pressure:" << PRESSURE << endl;

            if (PRESSURE < Pc && Rd < Rm)
            {

                // cout<<"ok1"<<endl;

                particle_id[i] = id_DV;

                grrateRd_mean = rv / (4 * M_PI * Rd * Rd);

                normal_distribution<double> dist7(grrateRd_mean, grrateRd_stdv);

                while (true)

                {
                    temp_Rd = dist7(gen);
                    if (temp_Rd > 0)
                        break;
                }

                // temp = normal2(grrateRd_mean, grrateRd_stdv);

                delta_Rd = temp_Rd * (double)time_step;

                // cout<<delta_Rd<<endl;

                radius_list[i] = Rd + delta_Rd;

                /*if (2 * Rd > temp_DIA)

                    temp_DIA = 2 * Rd;*/

                gammaMt_vsc[i] = eta * radius_list[i] * 6 * M_PI;
            }

            else if (PRESSURE < Pc && Rd >= Rm)
            {

                // cout << i << "Radius:" << Rd << "Pressure:" << PRESSURE << endl;

                // particle_id[i] = id_DV;

                numParticles++;

                unit_vector();

                // cout<<uv[0]<<" "<<uv[1]<<" "<<uv[2]<<endl;

                Rdtr = Rd * pow(2, -1.0 / 3.0);

                // cout<<Rd<<" "<<Rdtr<<endl;

                temp_DIST = Rd - Rdtr;

                // cout << "Rdc2" << endl;

                // cout << numParticles << endl;

                radius_list.push_back(0.0);

                gammaMt_vsc.push_back(0.0);

                nu_list.push_back(0.0);

                E_list.push_back(0.0);

                c1_list.push_back(0.0);

                c2_list.push_back(0.0);

                pos_x.push_back(0.0);

                pos_y.push_back(0.0);

                pos_z.push_back(0.0);

                particle_id.push_back(0.0);

                int counter1 = 0;

                for (int m = i; m < numParticles; m += (numParticles - i - 1))

                {

                    counter1++;

                    radius_list[m] = Rdtr;

                    gammaMt_vsc[m] = eta * Rdtr * 6 * M_PI;

                    nu_list[m] = dist2(gen);

                    while (true)

                    {

                        temp3 = dist3(gen);

                        if (temp3 > 0)

                            break;
                    }

                    E_list[m] = temp3;

                    while (true)

                    {

                        temp4 = dist4(gen);

                        if (temp4 >= 0 && temp4 <= 1)

                            break;
                    }

                    c1_list[m] = temp4;

                    while (true)

                    {

                        temp5 = dist4(gen);

                        if (temp5 >= 0 && temp5 <= 1)

                            break;
                    }

                    c2_list[m] = temp5;
                }

                /*if (counter1 != 2)

                    cout << "error" << endl;*/

                // cout<<temp<<endl;

                particle_id[i] = id_DV;

                particle_id[numParticles - 1] = id_lv;

                temp_x = pos_x[i];

                temp_y = pos_y[i];

                temp_z = pos_z[i];

                uniform_int_distribution<int> dist9(0, 1);

                temp = dist9(gen);

                // cout << temp << endl;

                if (temp == 0)
                {

                    pos_x[i] = temp_x + temp_DIST * uv[0];

                    // pos_x[i] = PBC(pos_x[i]);

                    pos_y[i] = temp_y + temp_DIST * uv[1];

                    // pos_y[i] = PBC(pos_[yi]);

                    pos_z[i] = temp_z + temp_DIST * uv[2];

                    // particle_id[i] = id_DV;

                    /*if(pos_x[numParticles-1]>eps && pos_x[numParticles-1]<(-eps))
                    cout<<"Error Pos_x"<<endl;*/

                    pos_x[numParticles - 1] = temp_x - temp_DIST * uv[0];

                    // pos_x[i] = PBC(pos_x[i]);

                    pos_y[numParticles - 1] = temp_y - temp_DIST * uv[1];

                    // pos_y[i] = PBC(pos_[yi]);

                    pos_z[numParticles - 1] = temp_z - temp_DIST * uv[2];

                    // particle_id[numParticles - 1] = id_lv;

                    /*double a1=pos_x[i]-pos_x[numParticles - 1];

                    double a2=pos_y[i]-pos_y[numParticles - 1];

                    double a3=pos_z[i]-pos_z[numParticles - 1];

                    temp_check=sqrt(a1*a1+a2*a2+a3*a3);

                    cout<<temp_check<<" "<<2*temp_DIST<<endl;*/
                }

                else

                {

                    pos_x[i] = temp_x - temp_DIST * uv[0];

                    // pos_x[i] = PBC(pos_x[i]);

                    pos_y[i] = temp_y - temp_DIST * uv[1];

                    // pos_y[i] = PBC(pos_[yi]);

                    pos_z[i] = temp_z - temp_DIST * uv[2];

                    /*if(pos_x[numParticles-1]>eps && pos_x[numParticles-1]<(-eps))
                    cout<<"Error Pos_x"<<endl;*/

                    pos_x[numParticles - 1] = temp_x + temp_DIST * uv[0];

                    // pos_x[i] = PBC(pos_x[i]);

                    pos_y[numParticles - 1] = temp_y + temp_DIST * uv[1];

                    // pos_y[i] = PBC(pos_[yi]);

                    pos_z[numParticles - 1] = temp_z + temp_DIST * uv[2];

                    /*double b1=pos_x[i]-pos_x[numParticles - 1];

                    double b2=pos_y[i]-pos_y[numParticles - 1];

                    double b3=pos_z[i]-pos_z[numParticles - 1];*/

                    // temp_check=sqrt(b1*b1+b2*b2+b3*b3);

                    // cout<<temp_check<<" "<<2*temp_DIST<<endl;
                }

                // cout << "Rdc5" << endl;
            }

            /*else
            {
                delta_Rd = 0;
                // cout<<"ok3"<<endl;
            }*/

            if (i == (N0_t - 1))
            {
                // cout << "rad function breaks at:" << i << endl;
                break;
            }
        }
    }
    /*{cout<<i<<endl;
        break;}*/
}

void time_advance()
{
    double temp, temp_DIST, temp_id;
    double temp_gamma, temp_gammaad1;
    double forceM;
    double dtbygammaMt, gammaMt;
    double prob_check;
    double Rds;

    // double grrateRd_mean;
    int counter2 = 0;
    int counter3 = 0;

    for (int i = 0; i < numParticles; i++)
    {

        temp_id = particle_id[i];

        if (temp_id == id_dd)
        // if (radius_list[i] < flag_2)
        {
            // C3.insert(i);
            continue;
        }

        else
        {
            // cout << prob_check << endl;

            // Rds=radius_list[i];

            prob_check = dist12(gen);

            if (prob_check < kdt)

            {

                // cout << "prob_check:" << prob_check << endl;

                // radius_list[i] = flag_1;

                particle_id[i] = id_dd;

                // cout << "id of dd particle at time advance function:" << i << endl;

                // C3.insert(i);

                // cout << radius_list[i] << endl;

                dd_numParticles++;

                continue;

                // cout << "continue does not work" << endl;
            }

            /*if (radius_list[i] < flag_2)

                continue;*/

            // rad_change(pressure[i],radius_list[i]);

            else

            {
                temp_gammaad1 = force_x[i] * gammaMt_adx1[i] + force_y[i] * gammaMt_ady1[i] + force_z[i] * gammaMt_adz1[i];

                // cout<<"ok"<<temp_gammaad1<<endl;

                // forceM = sqrt(force_x[i] * force_x[i] + force_y[i] * force_y[i] + force_z[i] * force_z[i]);

                /*if (forceM < eps)
                {
                    // cout << forceM << endl;
                    counter3++;
                }*/

                if (gammaMt_ad2[i] < eps)
                {
                    gammaMt = gammaMt_vsc[i];
                    // counter2++;
                }

                else

                {
                    if (temp_gammaad1 < eps && temp_gammaad1 > (-eps))

                        gammaMt = gammaMt_vsc[i] + 0.25 * gamma_max * gammaMt_ad2[i];

                    else

                    {
                        forceM = sqrt(force_x[i] * force_x[i] + force_y[i] * force_y[i] + force_z[i] * force_z[i]);

                        if (forceM < eps)

                            cout << "error" << endl;

                        temp_gammaad1 = temp_gammaad1 / forceM;

                        gammaMt = gammaMt_vsc[i] + 0.25 * gamma_max * (gammaMt_ad2[i] + temp_gammaad1);
                    }
                }

                dtbygammaMt = time_step / (gammaMt);

                dx = dtbygammaMt * force_x[i];

                pos_x[i] = pos_x[i] + dx;

                /*if (temp >= boxDimension)
                    pos_x[i] = temp - boxDimension;
                else if (temp < 0.0)
                    pos_x[i] = temp + boxDimension;
                else
                    pos_x[i] = temp;*/

                // pos_x[i] = PBC(temp);

                /*if(forceM==0.0)

                {counter2++;
                cout<<temp_gammaad1<<" "<<gammaMt<<" "<<dtbygammaMt<<" "<<temp<<" "<<pos_x[i]<<endl;}*/

                dy = dtbygammaMt * force_y[i];

                pos_y[i] = pos_y[i] + dy;

                /*if (temp >= boxDimension)
                    pos_y[i] = temp - boxDimension;
                else if (temp < 0.0)
                    pos_y[i] = temp + boxDimension;
                else
                    pos_y[i] = temp;*/

                // pos_y[i] = PBC(temp);

                dz = dtbygammaMt * force_z[i];

                pos_z[i] = pos_z[i] + dz;

                /*if (temp >= boxDimension)
                    pos_z[i] = temp - boxDimension;
                else if (temp < 0.0)
                    pos_z[i] = temp + boxDimension;
                else
                    pos_z[i] = temp;*/

                // pos_z[i] = PBC(temp);

                // temp = (dx * dx + dy * dy + dz * dz);

                // cout << "position change done" << endl;

                // cout << "force calc done" << endl;

                /*if (temp > maxD)
                    maxD = temp;*/
            }
        }
    }

    radius_change();

    // cout << "radius change done" << endl;

    calculate_force();
    /*nebListCounter++;
    if (2.0 * ((double)nebListCounter) * sqrt(maxD) > DR * MAXDIA)
        updateNebzLists();*/
    /*if (counter2 != counter3)
        cout << "Error" << endl;*/

    // return;
}

/*
#######

void datasavingsvectorDynamic(int vectorSize, vector<double> &vec1, vector<double> &vec2)
{

    vec1.assign(vectorSize, 0.0);
    vec2.assign(vectorSize, 0.0);

    for (int i = 0; i < vectorSize; i++)

        vec1[i] = (double)(i * time_step);

    // cout << vec1[0] << endl;

    // cout << "ok" << endl;
}*/

/*
############

void datatofile(string &identity, vector<double> &data1, vector<double> &data2)
{

    string filename1 = identity + "_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(ensemblesize) + "_ens_" + to_string(filenum) + "_filenum.txt";



    ofstream fout;

    fout.open(filename1.c_str());

    for (int i = 0; i < data1.size(); i++)

    {

        fout << data1[i] << " " << data2[i] / (double)ensemblesize << endl;
    }

    fout.close();
}

*/

void data_container_clear_1()
{
    /*vector<double>().swap(pos_x0);
    vector<double>().swap(pos_y0);
    vector<double>().swap(pos_z0);*/

    vector<double>().swap(uv);
    vector<double>().swap(pos_x);
    vector<double>().swap(pos_y);
    vector<double>().swap(pos_z);
    vector<double>().swap(force_x);
    vector<double>().swap(force_y);
    vector<double>().swap(force_z);
    vector<double>().swap(pressure);
    vector<double>().swap(radius_list);
    vector<double>().swap(E_list);
    vector<double>().swap(gammaMt_vsc);
    vector<double>().swap(nu_list);
    vector<double>().swap(c1_list);
    vector<double>().swap(c2_list);
    vector<double>().swap(gammaMt_adx1);
    vector<double>().swap(gammaMt_ady1);
    vector<double>().swap(gammaMt_adz1);
    vector<double>().swap(gammaMt_ad2);
    vector<int>().swap(particle_id);
}

void dynamicdata_container_clear_1()
{

    vector<double>().swap(TIME);
    // vector<int>().swap(timePRINT_indx);
    // set<int>().swap(timePRINT_indxSet);
    vector<int>().swap(timePRINT_chck);
}

void dynamicdata_container_clear_2()
{

    vector<int>().swap(PARTICLE_NUM);
    // vector<double>().swap(MSD);
    //  vector<double>().swap(CORR_R);
    //  vector<double>().swap(CORR_FUNC);
}

/*void data_container_clear_2()
{

    vector<double>().swap(MSD_time);
    vector<double>().swap(TIME);
    vector<double>().swap(CORR_func);
    vector<double>().swap(CORR_R);
}*/

int main(int argc, char *argv[])
{
    // boxDimension = 10.0;

    auto start = high_resolution_clock::now();

    N0 = atoi(argv[1]);
    fad = atof(argv[2]);
    //boxDimension = atof(argv[3]);
    stopTime = atoi(argv[3]);
    ensemblesize = atoi(argv[4]);
    filenum = atoi(argv[5]);

    // cout<<stopTime<<" "<<time_step<<endl;
    DynamicVectorDim = ((stopTime - startTime) / time_step) + 1;
    cout << "DynamicVectorDim:" << DynamicVectorDim << endl;

    cout << "L:" << boxDimension << " "
    << "N0:" << N0 << " "
    << "rv:" << rv << " "
    << "Kdt:" << kdt << " " << endl;
    //<< "flag_1:" << flag_1 << " "
    //<< "flag_2:" << flag_2 << endl;

    // int vector1dim = ((stopTime - startTime) / time_step) + 1;

    // cout << "vector1dim:"<<vector1dim << endl;

    // datasavingsvectorDynamic(vector1dim, TIME, PARTICLE_NUM);

    // cout << "vector1dim:" << TIME.size() << endl;

    DynamicDataSavings1();

    // cout<<"DDcontianerok"<<endl;

    // cout<<TIME.size()<<endl;

    for (int l = 0; l < ensemblesize; l++)
    {
        // data_container_initialize_1();

        int p_count = 0;

        //double

        initialize_system();

        // cout << "ok" << endl;
        cout << "Initial Phi: " << phi << endl;

        DynamicDataSavings2();

        // inititalDataPrint();

        ofstream fout1, fout2,fout3,fout4;

        //string filename1 = "lv_Position_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.txt";

        string filename1 = "lv_Position_" + to_string(N0) + "_N0_"+ to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.dat"; 

        //string filename2 = "dd_Position_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.txt";

        string filename2 = "dd_Position_" + to_string(N0) + "_N0_"+ to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.dat";    

        //string filename2 = "Particle_Num_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.txt";

       // string filename3 = "lv_Cell_Properties_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.txt";

        string filename3 = "lv_Cell_Properties_" + to_string(N0) + "_N0_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.dat";

        //string filename4 = "dd_Cell_Properties_" + to_string(N0) + "_N0_" + to_string(boxDimension) + "_L_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.txt";

        string filename4 = "dd_Cell_Properties_" + to_string(N0) + "_N0_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(l + 1) + "_ens_" + to_string(filenum) + "_filenum.dat";

        fout1.open(filename1.c_str());

        fout2.open(filename2.c_str());

        fout3.open(filename3.c_str());

        fout4.open(filename4.c_str());

        for (int i = 0; i < TIME.size(); i++)
        {

            PARTICLE_NUM[i] = numParticles - dd_numParticles;

            // Condition (a): Print data for i < 10
            if (i < 10)
            {
                
                fout1<<"ITEM: TIMESTEP\n";
                fout1<<i<<"\n";
                fout1<<"ITEM: NUMBER OF ATOMS\n";
                fout1<<PARTICLE_NUM[i]<<"\n";
                fout1<<"ITEM: BOX BOUNDS  pp pp pp\n";
                for(int z=0;z<3;z++)
                fout1<<-200.0<<'\t'<<200.0<<"\n";
                //fout1<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
                fout1<<"ITEM: ATOMS x y z radius type\n";

                fout3<<"ITEM: TIMESTEP\n";
                fout3<<i<<"\n";
                fout3<<"ITEM: ATOMS nu E c_1 c_2 type\n";

                fout2<<"ITEM: TIMESTEP\n";
                fout2<<i<<"\n";
                fout2<<"ITEM: NUMBER OF ATOMS\n";
                fout2<<dd_numParticles<<"\n";
                fout2<<"ITEM: BOX BOUNDS  pp pp pp\n";
                for(int z=0;z<3;z++)
                fout2<<-200.0<<'\t'<<200.0<<"\n";
                //fout2<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
                fout2<<"ITEM: ATOMS x y z radius type\n";

                fout4<<"ITEM: TIMESTEP\n";
                fout4<<i<<"\n";
                fout4<<"ITEM: ATOMS nu E c_1 c_2 type\n";

                for (int m = 0; m < numParticles; m++)
                {
                    // fout1 << fixed << setprecision(10) << m << '\t' << pos_x[m] << '\t' << pos_y[m] << '\t' << pos_z[m] << '\t' << radius_list[m] << '\t' << particle_id[m] << endl;
                    if(particle_id[m]!=id_dd)

                    {

                        /*fout1<<0<<'\t'<<129<<endl;
                        fout1<<0<<'\t'<<129<<endl;*/
                        fout1 << fixed << setprecision(10) << pos_x[m] << '\t'
                              << pos_y[m] << '\t'
                              << pos_z[m] << '\t'
                              << radius_list[m] << '\t'
                              << m << endl;

                        fout3 << fixed << setprecision(10) << nu_list[m] << '\t'
                              << E_list[m] << '\t'
                              << c1_list[m] << '\t'
                              << c2_list[m] << '\t'
                              << m << endl;

                    }          

                    else
                    {

                        //fout2<<"ITEM:"<< '\t'<<"TIMESTEP"<<endl;
                        //fout2<<fixed<<setprecision(10)<<TIME[i]<<'\t'<<dd_numParticles<<endl;
                        
                        fout2 << fixed << setprecision(10) << pos_x[m] << '\t'
                              << pos_y[m] << '\t'
                              << pos_z[m] << '\t'
                              << radius_list[m] << '\t'
                              << m << endl;

                        fout4 << fixed << setprecision(10) << nu_list[m] << '\t'
                              << E_list[m] << '\t'
                              << c1_list[m] << '\t'
                              << c2_list[m] << '\t'
                              << m << endl;

                    }     

                }
                // PARTICLE_NUM[i] = (double)(numParticles - dd_numParticles);
                //fout2 << fixed << setprecision(10) << TIME[i] << '\t' << PARTICLE_NUM[i] << endl;

                cout << "time step: " << i << " is complete" << endl;
                time_advance();
            }
        
        else if (i >= 10 && i < 7e5)
        {

            p_count++;

            if (timePRINT_chck[i] == 1)
            {
                
                fout1<<"ITEM: TIMESTEP\n";
                fout1<<i<<"\n";
                fout1<<"ITEM: NUMBER OF ATOMS\n";
                fout1<<PARTICLE_NUM[i]<<"\n";
                fout1<<"ITEM: BOX BOUNDS  pp pp pp\n";
                for(int z=0;z<3;z++)
                fout1<<-200.0<<'\t'<<200.0<<"\n";
                //fout1<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
                fout1<<"ITEM: ATOMS x y z radius type\n";

                fout3<<"ITEM: TIMESTEP\n";
                fout3<<i<<"\n";
                fout3<<"ITEM: ATOMS nu E c_1 c_2 type\n";

                fout2<<"ITEM: TIMESTEP\n";
                fout2<<i<<"\n";
                fout2<<"ITEM: NUMBER OF ATOMS\n";
                fout2<<dd_numParticles<<"\n";
                fout2<<"ITEM: BOX BOUNDS  pp pp pp\n";
                for(int z=0;z<3;z++)
                fout2<<-200.0<<'\t'<<200.0<<"\n";
                //fout2<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
                fout2<<"ITEM: ATOMS x y z radius type\n";

                fout4<<"ITEM: TIMESTEP\n";
                fout4<<i<<"\n";
                fout4<<"ITEM: ATOMS nu E c_1 c_2 type\n";

                for (int m = 0; m < numParticles; m++)
                {
                        // fout1 << fixed << setprecision(10) << m << '\t' << pos_x[m] << '\t' << pos_y[m] << '\t' << pos_z[m] << '\t' << radius_list[m] << '\t' << particle_id[m] << endl;
                    if(particle_id[m]!=id_dd)

                    {

                        fout1 << fixed << setprecision(10) << pos_x[m] << '\t'
                              << pos_y[m] << '\t'
                              << pos_z[m] << '\t'
                              << radius_list[m] << '\t'
                              << m << endl;

                        fout3 << fixed << setprecision(10) << nu_list[m] << '\t'
                              << E_list[m] << '\t'
                              << c1_list[m] << '\t'
                              << c2_list[m] << '\t'
                              << m << endl;
                
                    }

                    else
                    {

                        fout2 << fixed << setprecision(10) << pos_x[m] << '\t'
                              << pos_y[m] << '\t'
                              << pos_z[m] << '\t'
                              << radius_list[m] << '\t'
                              << m << endl;

                        fout4 << fixed << setprecision(10) << nu_list[m] << '\t'
                              << E_list[m] << '\t'
                              << c1_list[m] << '\t'
                              << c2_list[m] << '\t'
                              << m << endl;

                    }        


                }
                    // PARTICLE_NUM[i] = (double)(numParticles - dd_numParticles);
                    //fout2 << fixed << setprecision(10) << TIME[i] << '\t' << PARTICLE_NUM[i] << endl;

                if (p_count % 5 == 0)

                    cout << "time step: " << i << " is complete"
                << "& Particle No. is" << PARTICLE_NUM[i] << endl;
            }

            time_advance();
        }

        else
            //(i == 1e5)
        {

            fout1<<"ITEM: TIMESTEP\n";
            fout1<<i<<"\n";
            fout1<<"ITEM: NUMBER OF ATOMS\n";
            fout1<<PARTICLE_NUM[i]<<"\n";
            fout1<<"ITEM: BOX BOUNDS  pp pp pp\n";
            for(int z=0;z<3;z++)
            fout1<<-200.0<<'\t'<<200.0<<"\n";
                //fout1<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
            fout1<<"ITEM: ATOMS x y z radius type\n";

            fout3<<"ITEM: TIMESTEP\n";
            fout3<<i<<"\n";
            fout3<<"ITEM: ATOMS nu E c_1 c_2 type\n";

            fout2<<"ITEM: TIMESTEP\n";
            fout2<<i<<"\n";
            fout2<<"ITEM: NUMBER OF ATOMS\n";
            fout2<<dd_numParticles<<"\n";
            fout2<<"ITEM: BOX BOUNDS  pp pp pp\n";
            for(int z=0;z<3;z++)
            fout2<<-200.0<<'\t'<<200.0<<"\n";
                //fout2<<"ITEM:"<<'\t'<<"ATOMS"<<'\t'<<"x"<<'\t'<<"y"<<'\t'<<"z"<<'\t'<<"radius"<<'\t'<<"type"<<'\t'<<"id"<<endl;
            fout2<<"ITEM: ATOMS x y z radius type\n";

            fout4<<"ITEM: TIMESTEP\n";
            fout4<<i<<"\n";
            fout4<<"ITEM: ATOMS nu E c_1 c_2 type\n";
            
           

            for (int m = 0; m < numParticles; m++)
            {
                // fout1 << fixed << setprecision(10) << m << '\t' << pos_x[m] << '\t' << pos_y[m] << '\t' << pos_z[m] << '\t' << radius_list[m] << '\t' << particle_id[m] << endl;
                
                if(particle_id[m]!=id_dd)

                {

                    fout1 << fixed << setprecision(10) << pos_x[m] << '\t'
                          << pos_y[m] << '\t'
                          << pos_z[m] << '\t'
                          << radius_list[m] << '\t'
                          << m << endl;

                    fout3 << fixed << setprecision(10) << nu_list[m] << '\t'
                          << E_list[m] << '\t'
                          << c1_list[m] << '\t'
                          << c2_list[m] << '\t'
                          << m << endl;
                
                }

                else
                {

                    fout2 << fixed << setprecision(10) << pos_x[m] << '\t'
                          << pos_y[m] << '\t'
                          << pos_z[m] << '\t'
                          << radius_list[m] << '\t'
                          << m << endl;

                    fout4 << fixed << setprecision(10) << nu_list[m] << '\t'
                          << E_list[m] << '\t'
                          << c1_list[m] << '\t'
                          << c2_list[m] << '\t'
                          << m << endl;

                }  
                                
            }
                // PARTICLE_NUM[i] = (double)(numParticles - dd_numParticles);
                //fout2 << fixed << setprecision(10) << TIME[i] << '\t' << PARTICLE_NUM[i] << endl;
            cout << "time step: " << i << " is complete" << endl;
                // time_advance();
        }

            // cout << "time step: " << i << " is complete" << endl;
            // time_advance();
    }

        /*for (int m = 0; m < numParticles; m++)

            fout << fixed << setprecision(10) << m << '\t' << pos_x[m] << '\t' << pos_y[m] << '\t' << pos_z[m] << '\t' << particle_id[m] << endl;

        PARTICLE_NUM[i] = (double)(numParticles - dd_numParticles);

        fout << fixed << setprecision(10) << TIME[i] << '\t' << PARTICLE_NUM[i] << endl;*/

        // if(i==1)
        // SnapshotDataPrint();

        /*cout << "step no.:" << i << " "
             << "No. of particles:" << numParticles << " "
             << "No. of dead particles: " << dd_numParticles << endl;*/

        /*cout << "time step: " << i << "is complete" << endl;

        time_advance();*/

        // PositionDataPrint();

    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();

    data_container_clear_1();

        // cout<<"inverse number density:"<<inv_num_density<<endl;

    cout << l << " no. experiment is done" << endl;

    dynamicdata_container_clear_2();

        // calc_dynamics();
}

dynamicdata_container_clear_1();

    /*

    ######

        string filetag_1 = "Particle Number";

        datatofile(filetag_1, TIME, PARTICLE_NUM);

    */

    // ofstream fout;

    /*fout.open("particles_position.txt");

    /
    for (int i = 0; i < numParticles; i++)
    {

        fout << i << " " << x_pos[i] << " " << y_pos[i] << endl;
    }*/

    //}

    // fout.close();

auto stop = high_resolution_clock::now();

auto duration = duration_cast<seconds>(stop - start);

string filename8 = "#Time_Record_for_" + to_string(N0) + "_N0_" + to_string(fad) + "_fad_" + to_string(stopTime) + "_time_" + to_string(ensemblesize) + "_ens_" + to_string(filenum) + "_filenum.dat";

ofstream fout8;

fout8.open(filename8.c_str());

fout8 << "#Time taken for this data file:" << duration.count() << "seconds" << endl;

fout8.close();
    // Rest of your code

    /*cout << "Does indexing work?" << endl;

    for (auto i : C1)
    {
        cout << i << endl;
    }
    cout << "C1" << endl;
    for (auto i : C2)
    {
        cout << i << endl;
    }
    cout << "C2" << endl;
    for (auto i : C3)
    {
        cout << i << endl;
    }
    cout << "C3" << endl;*/
return 0;
}
