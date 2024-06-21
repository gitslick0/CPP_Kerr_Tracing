# include <iostream>
# include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>

extern "C" {
    void __blcoordinate_MOD_center_of_image(const double* robs, const double* sinobs, const double* muobs, const double* a_spin, const double* scal,
                                            const double* velocity, double* alphac, double* betac);
    void __blcoordinate_MOD_lambdaq(const double* alpha, const double* beta, const double* robs, const double* sinobs, const double* muobs, const double* a_spin, const double* scal,
                                    const double* velocity, double* f1234, double* lambda, double* q);
    double __blcoordinate_MOD_rms(const double* a_spin);
    void __blcoordinate_MOD_radiustp(const double* f12341,const double* a_spin,const double* robs, const double* lambda, const double* q,
                                    double* r_tp1, double* r_tp2, int* reals, bool* robs_eq_rtp, bool* indrhorizon, int* cases, std::complex<double>* bb);
    double __blcoordinate_MOD_r2p(const double* f1234r, const double* rend, const double* lambda, const double* q, const double* a_spin, const double* robs,
                                const double* scal, const int* t1, const int* t2);
 void __blcoordinate_MOD_ynogk(const double* p, const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin,
                                const double* robs, const double* scal, const double* radi, const double* mu, const double* phi, const double* time, const double* sigma);
}

// Convert external functions
std::array<double, 2> fcenter_of_image(const double robs, const double theta, const double a_spin, std::array<double, 3> velocity){
    double alphac = 0.0, betac = 0.0;
    double sinobs = std::sin(theta), muobs = std::cos(theta);
    double scal = 1.0;
    __blcoordinate_MOD_center_of_image(&robs, &sinobs, &muobs, &a_spin, &scal, velocity.data(), &alphac, &betac);
    std::array<double,2> rvec = {static_cast<float>(alphac), static_cast<float>(betac)};
    return rvec;
}

std::tuple<std::array<double, 4>, double, double> flambdaq(const std::array<double, 2> alphabeta, const double robs, const double theta, const double a_spin,
const std::array<double, 3> velocity){
    std::array<double, 4> f1234;
    double lambda, q;
    double sinobs = std::sin(theta), muobs = std::cos(theta), scal = 1.0;
    __blcoordinate_MOD_lambdaq(&alphabeta.data()[0], &alphabeta.data()[1], &robs, &sinobs, &muobs, &a_spin, &scal, velocity.data(), f1234.data(), &lambda, &q);

    return std::make_tuple(f1234, lambda, q);
}

std::tuple<double, double, int, bool, bool, int, std::array<std::complex<double>, 4>> fradiustp(const std::array<double, 4> f1234, const double a_spin, const double robs,
                                                                                const double lambda, const double q){
    double r_tp1, r_tp2;
    int reals, cases;
    bool robs_eq_rtp, indrhorizon;
    std::array<std::complex<double>, 4> bb;

    __blcoordinate_MOD_radiustp(f1234.data(), &a_spin, &robs, &lambda, &q, &r_tp1, &r_tp2, &reals, &robs_eq_rtp, &indrhorizon, &cases, bb.data());
    return std::make_tuple(r_tp1, r_tp2, reals, robs_eq_rtp, indrhorizon, cases, bb);
}

double frms(const double a_spin){
    return __blcoordinate_MOD_rms(&a_spin);
}

double fr2p(std::array<double, 4> f1234, const double rend, const double lambda, const double q, const double a_spin, const double robs,
             const int t1, const int t2){
    double scal = 1.0, f1234r = f1234.data()[0];
    return __blcoordinate_MOD_r2p(&f1234r, &rend, &lambda, &q, &a_spin, &robs, &scal, &t1, &t2);
}

std::array<double, 5> YNOGK(double p, std::array<double, 4> f1234, std::array<double, 2> Photon_motion_constants,std::array<double, 3> source_pos, double a_spin){
  double scal = 1.0;
  double radi = 0.0, mu = 0.0, phi = 0.0, time = 0.0, sigma = 0.0;
  __blcoordinate_MOD_ynogk(&p, f1234.data(), &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &source_pos.data()[1],
                                   &source_pos.data()[2], &a_spin, &source_pos.data()[0], &scal, &radi, &mu, &phi, &time, &sigma);
  return {radi, mu, phi, time, sigma};
}



int main(){
    const int height = 1000, width =1000;
    std::array<double, 4> f1234;
    double lambda, q, sigmaa;
    std::tuple<double, double, int, bool, bool, int, std::array<std::complex<double>, 4>> frtp;
    std::array<double, 5> ynogkres;
    Eigen::MatrixXf diskr(height, width);
    //std::vector<float> data_array(width*height);
    //std::vector<float> alphas(width*height);
    //std::vector<float> betas(width*height);
    // Example input data
    float a_spin = 0.998, c_obs = M_PI/2, robs = float(1E6), scal = 1.0;
    std::array<double, 3> velocity = {0.0, 0.0, 0.0};
    std::array<double, 2> coi = {1.0, 1.0};
    // Size of Image
    float deltax = 12.0/width, deltay = 12.0/height;
    // Set constants
    float x = 43.0, y = float(1E70);
    int t1 = 0, t2 = 0;
    float p = 1.0, rhorizon = 1+sqrt(1-a_spin+a_spin);

    // Find center of image
    coi = fcenter_of_image(robs, c_obs, a_spin, velocity);
    std::cout << "center of image ";
    for(auto alp : coi){
        std::cout << alp << " ";
    }
    std::cout << std::endl;


    // Loop over pixels in image
    for(int ii = 0; ii < height; ii++){
        sigmaa = 0.0;
        float beta = 6.0 - ii*deltay + coi[1];
        for(int jj = 0; jj < width; jj++){
            float alpha = 4.0 - jj*deltax + coi[0];
            std::tuple<std::array<double, 4>, double, double> f1234lambdaq = flambdaq({alpha, beta}, robs, c_obs, a_spin, velocity);
            f1234 = std::get<0>(f1234lambdaq);
            lambda = std::get<1>(f1234lambdaq);
            q = std::get<2>(f1234lambdaq);
            frtp = fradiustp(f1234, a_spin, robs, lambda, q);
            double r_tp1 = std::get<0>(frtp), r_tp2 = std::get<1>(frtp);
            if(r_tp1 <= rhorizon){
                t1 = 0, t2 = 0;
                p = fr2p(f1234, rhorizon, lambda, q, a_spin, robs, t1, t1);
                ynogkres = YNOGK(p,f1234,{lambda,q},{robs, std::sin(c_obs),std::cos(c_obs)},a_spin);
                sigmaa = static_cast<float>(ynogkres[4]);
                //data_array[ii*width + jj] =  isnan(sigmaa) ? 0.0 : static_cast<float>(sigmaa);
                diskr(jj, ii) = isnan(sigmaa) ? 0.0 : std::abs(sigmaa);
            }
            else{
                t1 = 1, t2 = 0;
                p = fr2p(f1234, rhorizon, lambda, q, a_spin, robs, t1, t1);
                ynogkres = YNOGK(p,f1234,{lambda,q},{robs, std::sin(c_obs),std::cos(c_obs)},a_spin);
                sigmaa = static_cast<float>(ynogkres[4]/2.0);
                //data_array[ii * width + jj] =  isnan(sigmaa) ? 0.0 : static_cast<float>(sigmaa);
                diskr(jj, ii) = isnan(sigmaa) ? 0.0 : std::abs(sigmaa);
            }
        //alphas[ii*width + jj] = alpha;
        //betas[ii*width + jj] = beta;
        }
    }

float maxr = diskr.maxCoeff();
float minr = diskr.minCoeff();
std::cout << "Max: " << maxr << ", Min: " << minr << std::endl;

// Convert Matrix to OpenCV
cv::Mat image(height, width, CV_32F, diskr.data());
cv::Mat displayImage;
cv::normalize(image, displayImage, 0, 255, cv::NORM_MINMAX, CV_8U);

// Display the Image
cv::namedWindow("Image", cv::WINDOW_AUTOSIZE);
cv::imshow("Image", displayImage);
cv::waitKey(0);




return 0;

}