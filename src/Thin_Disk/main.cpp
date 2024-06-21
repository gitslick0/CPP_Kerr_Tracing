# include <iostream>
# include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>
#include <algorithm>

extern "C" {
    void __blcoordinate_MOD_metricg(const double* robs, const double* sinobs, const double* muobs, const double* a_spin,
                                    double* somiga, double* expnu, double* exppsi, double* expmu1, double* expmu2);
    double __blcoordinate_MOD_rms(const double* a_spin);
    void __blcoordinate_MOD_center_of_image(const double* robs, const double* sinobs, const double* muobs, const double* a_spin, const double* scal, const double* velocity,
                                            double* alpha_c, double* beta_c);
    void __blcoordinate_MOD_lambdaq(const double* alpha, const double* beta, const double* robs, const double* sinobs, const double* muobs, const double* a_spin,
                                    const double* scal, const double* velocity, double* f1234, double* lambda, double* q);
    double __blcoordinate_MOD_pemdisk_all(const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs,
                                        const double* a_spin, const double* robs, const double* scal, const double* mu, const double* rout, const double* rin);
    double __blcoordinate_MOD_radius(const double* p, const double* f1234r, const double* lambda, const double* q, const double* a_spin, const double* robs, const double* scal);

}

int main(){

// Example Input Data:
const double dtor = M_PI/180.0;
const double a_spin = 0.998, theta = 90.0*dtor, sinobs = std::sin(theta), muobs = std::cos(theta), robs = 100; // theta = angle of observer, robs = radius of observer
const double scal = 1.0, theta_disk = 0.0, rdisk_out = 80.0, mudisk = std::cos((90.0 - theta_disk))*dtor; // theta_disk = inclination angle of disk, rdisk_out = outer edge radius of disk
int m = 1200; // image size

// Variables needed:
double somiga_obs, expnu_obs, exppsi_obs, expmu1_obs, expmu2_obs, rms1;
double somiga_em, expnu_em, exppsi_em, expmu1_em, expmu2_em;
std::array<double, 3> velocity;
std::array<double,4> f1234;
double alpha_c, beta_c, alpha, beta;
double delta_x, delta_y;
double lambda, q;
double pem, re;
double one = 1.0, zero = 0.0;
double bomiga, ut_em, g;
std::vector<float> g_vector(m*m);


__blcoordinate_MOD_metricg(&robs, &sinobs, &muobs, &a_spin, &somiga_obs, &expnu_obs, &exppsi_obs, &expmu1_obs, &expmu2_obs);
rms1 = __blcoordinate_MOD_rms(&a_spin);
velocity[0] = robs/(pow(robs, 1.5) + a_spin) * 0.0;
velocity[1] = robs/(pow(robs, 1.5) + a_spin) * std::cos(M_PI/4)*0.0;
velocity[2] = robs/(pow(robs, 1.5) + a_spin) * std::sin(M_PI/4)*0.0;

// calculate center of image
__blcoordinate_MOD_center_of_image(&robs, &sinobs, &muobs, &a_spin, &scal, velocity.data(), &alpha_c, &beta_c);
delta_x = 60.0/m;
delta_y = 28.0/m;

for(int ii = 0; ii < m; ++ii){
    beta=beta_c-ii*delta_y+13.0;
    for(int jj = 0; jj < m; ++jj){
        alpha=alpha_c-jj*delta_x+30.0;
        __blcoordinate_MOD_lambdaq(&alpha, &beta, &robs, &sinobs, &muobs, &a_spin, &scal, velocity.data(), f1234.data(), &lambda, &q);
        pem = __blcoordinate_MOD_pemdisk_all(f1234.data(), &lambda, &q, &sinobs, &muobs, &a_spin, &robs, &scal,
                                            &mudisk, &rdisk_out, &rms1);
        if(pem != -1.0 && pem != -2.0){
            re = __blcoordinate_MOD_radius(&pem, &f1234.data()[0], &lambda, &q, &a_spin, &robs, &scal);
            bomiga=1.0/(a_spin+pow(re, 1.5));
            __blcoordinate_MOD_metricg(&re,&one,&zero,&a_spin,&somiga_em,&expnu_em,&exppsi_em,&expmu1_em,&expmu2_em);
            ut_em=1.0/expnu_em/std::sqrt(1.0-pow(exppsi_em/expnu_em*(bomiga-somiga_em), 2)); 
            g = (1.0-somiga_obs*lambda)/expnu_obs/f1234[3]/(1.0-bomiga*lambda)/ut_em;
        }
        else{
            if(pem != -1.0){
                g = 0.0;
            }
            else{
                g = 0.1;
            }
        }
        g_vector[ii*m + jj] = g;
    }
}

// Prepare data:

auto minmax = std::minmax_element(g_vector.begin(), g_vector.end());
float minr = *minmax.first;
float maxr = *minmax.second;

// Print minimum and maximum value
std::cout << "Minimum " << minr << " Maximum " << maxr << std::endl;

for(auto &val : g_vector){
    if(val == 0.0){
        val = maxr;
    }
    else if(val == 0.1){
        val = maxr;
    }
    val = std::abs(val);
}

// Normalize and scale data:
std::vector<uint8_t> byteDiskr(m*m);
for(size_t ii = 0; ii < g_vector.size(); ++ii){
    byteDiskr[ii] = static_cast<uint8_t>((g_vector[ii] - minr) / (maxr - minr) * 255);
}

// Create OpenCV matrix
cv::Mat image(m, m, CV_8UC1, byteDiskr.data());
cv::applyColorMap(image, image, cv::COLORMAP_TWILIGHT_SHIFTED);
cv::flip(image, image, 0);

// Display Image
cv::imshow("Disk Image", image);
cv::waitKey(0);


    return 0;
}