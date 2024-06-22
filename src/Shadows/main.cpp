#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <complex>
#include <algorithm>
#include <opencv2/opencv.hpp>

extern "C" {
    void __blcoordinate_MOD_center_of_image(const double* robs, const double* sinobs, const double* muobs, const double* a_spin, const double* scal, const double* velocity,
                                            double* alpha_c, double* beta_c);
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

int main(){
    // Example input data
    double a_spin = 0.998, dtor = M_PI/180.0, theta = 90.0*dtor, robs = 1.0*pow(10,6), scal = 1.0;
    double sinobs = std::sin(theta), muobs = std::cos(theta);
    std::array<double, 3> velocity = {0.0, 0.0, 0.0};

    // Data we need later
    double alpha_c, beta_c, alpha, beta;
    double lambda, q;
    std::array<double, 4> f1234;
    int m = 1000; // size of image in pixels
    // results of radiustp
    double rtp1, rtp2;
    int reals;
    bool robs_eq_rtp, indrhorizon;
    int cases_of_tp;
    std::array<std::complex<double>, 4> bb;
    // results of ynogk
    double ra, mua, phia, timea, sigmaa;
    // data vector
    std::vector<float> sigmaa_vector(m*m, 0.0f);

    // set some constants initially
    int t1 = 0, t2;
    double p = 1.0, x = 43.0, y = 1.0E70;
    double rhorizon = 1.0 + std::sqrt(1.0-pow(a_spin,2));
    std::cout << " rhorizon = " << rhorizon << std::endl;

    // Calculate center of image (centered on black hole)
    std::cout << "Calculating center of image..." << std::endl;
    __blcoordinate_MOD_center_of_image(&robs, &sinobs, &muobs, &a_spin, &scal, velocity.data(), &alpha_c, &beta_c);
    double delta_x = 12.0/m;
    double delta_y = 12.0/m;
    std::cout << "Center of image calculated: alpha_c = " << alpha_c << ", beta_c = " << beta_c << std::endl;

    // Loop over pixels to calculate sigmaa
    for(int ii = 0; ii < m; ++ii){
        beta = 6.0 - ii*delta_y + beta_c;
        for(int jj = 0; jj < m; ++jj){
            alpha = 4.0 - jj*delta_x + alpha_c;
            __blcoordinate_MOD_lambdaq(&alpha,&beta,&robs,&sinobs,&muobs,&a_spin,&scal,velocity.data(),f1234.data(),&lambda,&q);
            __blcoordinate_MOD_radiustp(&f1234.data()[0],&a_spin,&robs,&lambda,&q,&rtp1,&rtp2,&reals,&robs_eq_rtp,&indrhorizon,&cases_of_tp,bb.data());

            if (ii * m + jj >= sigmaa_vector.size()) {
                std::cerr << "Index out of bounds: " << ii * m + jj << std::endl;
                return 1;
            }
            if(rtp1 <= rhorizon){
                t1 = 0, t2 = 0;
                p = __blcoordinate_MOD_r2p(f1234.data(),&rhorizon,&lambda,&q,&a_spin,&robs,&scal,&t1,&t2);
                __blcoordinate_MOD_ynogk(&p,f1234.data(),&lambda,&q,&sinobs,&muobs,&a_spin,&robs,&scal,&ra,&mua,&phia,&timea,&sigmaa);
                sigmaa_vector[ii*m + jj] = isnan(sigmaa) ? 0.90909 : sigmaa;
            } else {
                t1 = 1, t2 = 0;
                p = __blcoordinate_MOD_r2p(f1234.data(),&rhorizon,&lambda,&q,&a_spin,&robs,&scal,&t1,&t2);
                __blcoordinate_MOD_ynogk(&p,f1234.data(),&lambda,&q,&sinobs,&muobs,&a_spin,&robs,&scal,&ra,&mua,&phia,&timea,&sigmaa);
                sigmaa_vector[ii * m + jj] = isnan(sigmaa) ? 0.90909 : sigmaa/2;
            }
        }
    }

std::cout << "set sigmaa vector 0" << std::endl;

    auto minmax = std::minmax_element(sigmaa_vector.begin(), sigmaa_vector.end());
    float minr = *minmax.first;
    float maxr = *minmax.second;

//    for(auto &value : sigmaa_vector){
//        if(value == 0.90909){
//            value = minr;
//        }
//    }

    // Print minimum and maximum value
    std::cout << "Minimum " << minr << " Maximum " << maxr << std::endl;

    // Normalize and scale data:
    std::vector<uint8_t> byteDiskr(m*m);
    for(size_t ii = 0; ii < sigmaa_vector.size(); ++ii){
        byteDiskr[ii] = static_cast<uint8_t>((sigmaa_vector[ii] - minr) / (maxr - minr) * 255);
    }

// Create OpenCV matrix
cv::Mat image(m, m, CV_8UC1, byteDiskr.data());
cv::applyColorMap(image, image, cv::COLORMAP_BONE);
cv::flip(image, image, 0);

// Display Image
cv::imshow("Disk Image", image);
cv::waitKey(0);

    return 0;
}
