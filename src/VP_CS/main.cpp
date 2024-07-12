// Calculate the velocity of sources based on crosssectionally related accretion disk element
# include <iostream>
# include <cmath>
# include <VPCS/helper_functions.h>
# include <random>
# include <algorithm>
# include<valarray>
# include <array>



extern "C" {
    double Fp(const double* p, const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin,
              const double* robs, const double* scal, const int* t1, const int* t2, const double* paras){ //f1234(1:4), paras(1:10);
                 double radi, mu, phi, time, sigma;
                 __blcoordinate_MOD_ynogk(p, f1234, lambda, q, sinobs, muobs, a_spin, robs, scal, &radi, &mu, &phi, &time, &sigma);
                 return mu;
    }
};



int main(){
    // Initialize some variables for later
    std::array<double, 4> f1234;
    double lambda, q;
    std::array<double, 3> velocity = {0.0, 0.0, 0.0};

    // Set Black Hole Spin parameter
    double a_spin;
    std::cout << "Set a_spin: " << std::endl;
    std::cin >> a_spin;
    //Get Input Source position
    double robs, cobs, dtor = M_PI/180.0;

    std::cout << "Set source radius: " << std::endl;
    std::cin >> robs; std::cout << "Set source inclination: " << std::endl;
    std::cin >> cobs;

    double sinobs = std::sin(cobs*dtor), muobs = std::cos(cobs*dtor);

    // Set accretion disk to be sampled
    double rms = __blcoordinate_MOD_rms(&a_spin), rout = 1000.0;
    int npoints = 10; //How many elements we want to try
    
    // Create logarithmic spacing of the disk elements
    std::vector<double> r_elements = LogarithmicSample(rms, rout, npoints);
    // Tested so far, seems about right.

    // Generate Uniform Emission Angles
    std::vector<double> RAND1 = generateRandomNumbers(npoints, 0, 2*M_PI);
    std::vector<double> RAND2 = generateRandomNumbers(npoints, 0, 1);
    std::vector<double> PP_ARR(npoints), PR_ARR(npoints), PT_ARR(npoints);
    for(int jj = 0; jj < RAND1.size(); ++jj){
        PT_ARR[jj] = std::sin(RAND1[jj])*std::cos(RAND2[jj]);
        PR_ARR[jj] = std::sin(RAND1[jj])*std::sin(RAND2[jj]);
        PP_ARR[jj] = std::cos(RAND2[jj]);
    }
    // Loop over all test disk elements
    for(auto &r_em : r_elements){
        double sin_em = 1.0, mu_em = 0.0;
        // Loop over emitted photons
        for(int jj = 0; jj < PP_ARR.size(); ++jj){
            // Calculate initial momenta
            __blcoordinate_MOD_initialdirection(&PR_ARR[jj], &PT_ARR[jj], &PP_ARR[jj], &sin_em, &mu_em, &a_spin, &r_em,  velocity.data(),
                                                &lambda, &q, f1234.data());
            double p_zero = 0.1;
            int t1=0, t2=0;
            std::array<double, 10> paras; paras.fill(0.0);
            double scal = 1.0, rin = robs-0.3*robs, rout = robs + 0.3*robs, muup = muobs + 0.3*muobs, mudown = muobs + 0.3*muobs, phy1 = 2*M_PI, phy2 = 0.0;
            int caserange = 1.0;
            bool bisection = false;
            double result = __pemfinding_MOD_pemfind(f1234.data(), &lambda, &q, &sinobs, &muobs, &a_spin, &robs, &scal, &rin, &rout, &muup, &mudown,
                                    &phy1, &phy2, &caserange, &Fp, paras.data(), &bisection);
            std::cout << "jj = " << jj;
            std::cout << " and result = " << result << std::endl;
        }
    }

    return 0;

}