#ifndef VPCS_HELPERS
#define VPCS_HELPERS

# include<vector>
# include<cmath>
# include<iostream>
# include<random>

extern "C" {
    // functions to be imported from ynogk.f90
    double __blcoordinate_MOD_rms(const double* a_spin);
    void __blcoordinate_MOD_initialdirection(const double* pr, const double* ptheta, const double* pphi, const double* sinobs,
                                            const double* muobs, const double* a_spin, const double* robs, const double* velocity, double* lambda, double* q, double* f1234);
    double __pemfinding_MOD_pemfind(const double* f1234,const double* lambda,const double* q, const double* sinobs, const double* muobs, const double* a_spin,
                                    const double* robs, const double* scal, const double* rin, const double* rout, const double* muup, const double* mudown,
                                    const double* phy1, const double* phy2, const int* caserange, double (*Fp)(const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const double*),
                                    const double* paras, const bool bisection);
    void __blcoordinate_MOD_ynogk(const double* p, const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin,
                                const double* robs, const double* scal, const double* radi, const double* mu, const double* phi, const double* time, const double* sigma);

};

// logarithmic sampling
std::vector<double> LogarithmicSample(double min, double max, int npoints){
    std::vector<double> points;
    points.reserve(npoints);
    // Validate Input
    if (min <= 0 || max <= 0){
        throw std::invalid_argument("Range values must be greater than zero for logarithmic scale");
    }
    if(npoints <= 1){
        points.push_back(min);
        return points;
        //throw std::invalid_argument("Number of points must be greater than 1.");
    }

    // Compute the logarithmic spacing
    double log_min = std::log(min);
    double log_max = std::log(max);
    double log_step = (log_max - log_min)/(npoints - 1);

    for(int ii = 0; ii < npoints; ++ii){
        double log_value = log_min + ii*log_step;
        points.push_back(std::exp(log_value));
    }

    return points;
}

// generate random numbers
std::vector<double> generateRandomNumbers(int M, double min, double max){
    // create random number engine
    std::random_device rd;  // Seed the random number engine
    std::mt19937 gen(rd()); // Mersenne Twister engine

    // Define the distribution range
    std::uniform_real_distribution<> dis(min, max);

    // Generate M numbers
    std::vector<double> random_numbers;
    random_numbers.reserve(M);
    for(int ii = 0; ii < M; ++ii){
        random_numbers.push_back(dis(gen));
    }

    return random_numbers;

}

/* class Photon{
    private:
        std::array<double, 4> f1234,
}

class Emission_Source{
    private:
        double robs, cobs, a_spin
        int nPh // Number of Photons to emit

    public:
        Emission_Source(double inp_robs, double inp_cobs, double inp_a_spin) : robs(inp_robs), cobs(inp_cobs), a_spin(inp_a_spin), nPh(0);
        std::vector<Photon> Emission();

}
*/














#endif