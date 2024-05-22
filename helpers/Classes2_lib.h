#ifndef _Classes2_lib_
#define _Classes2_lib_

#include <iostream>
#include"coordinates.h"
#include"Drawables.h"
#include <cmath>
#include <ctime>
#include <iterator>
#include <array>
#include <vector>
#include <cstdlib>
#include <functional>
#include <GL/glut.h>
#include <GL/freeglut.h>

extern "C" {
  // functions I want to import
  double __blcoordinate_MOD_rms(const double *a_spin);
  void __blcoordinate_MOD_initialdirection(const double* pr, const double* ptheta, const double* pphi, const double* sinobs, const double* muobs, const double* a_spin, const double* r_obs,
                                           const double* velocity, const double* lambda, const double* q, const double* f1234);
  void __blcoordinate_MOD_ynogk(const double* p, const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin,
                                const double* robs, const double* scal, const double* radi, const double* mu, const double* phi, const double* time, const double* sigma);
  double __blcoordinate_MOD_p_total(const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin, const double* robs, const double* scal);
  double __blcoordinate_MOD_pemdisk(const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin, const double* robs,
                                    const double* scal, const double* mu, const double* rout, const double* rin);
  void __blcoordinate_MOD_metricg(const double* robs, const double* sinobs, const double* muobs, const double* a_spin, const double* somega, const double* expnu,
                                const double* exppsi, const double* expmu1, const double* expmu2);
  double __blcoordinate_MOD_mu2p(const double* f12343, const double* f12342, const double* lambda, const double* q, const double* mu, const double* sinobs, const double* muobs, const double* a_spin,
                                 const int* t1, const int* t2, const double* scal);
}

// Convert fortran wrappers to C++ functions

double frms(const double a_spin){
  return __blcoordinate_MOD_rms(&a_spin);
}
CoordVec4 finitial_direction (MomentumVector momenta, CoordVec3 source_pos, const double a_spin, CoordVec3 source_vel){
  double lambda = 0, q=0;
  std::array<double, 4> temp_f1234
  __blcoordinate_MOD_initialdirection(&momenta.data()[1], &momenta.data()[2], &momenta.data()[0], &source_pos.data()[1], &source_pos.data()[2], &a_spin, &source_pos.data()[0],
                                      source_vel.data(), &lambda, &q, tempf1234.data());
  return CoordVec4(f1234[0], f1234[1], f1234[2], f1234[3]);
}

std::array<double,2> calculate_photon_motion_constants (std::array<double, 3> momenta, std::array<double, 3> source_pos, const double a_spin, std::array<double, 3> source_vel){
  double lambda = 0; double q = 0;
  std::array<double,4> f1234;
  __blcoordinate_MOD_initialdirection(&momenta.data()[1], &momenta.data()[2], &momenta.data()[0], &source_pos.data()[1], &source_pos.data()[2], &a_spin, &source_pos.data()[0],
                                      source_vel.data(), &lambda, &q, f1234.data());
  return {lambda, q};
}


double fp_total(const std::array<double, 4> f1234, const std::array<double, 2> Photon_motion_constants, const std::array<double, 3> source_pos, const double a_spin){
  double scal = 1.0;
  return __blcoordinate_MOD_p_total(f1234.data(), &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &source_pos.data()[1], &source_pos.data()[2],
                                    &a_spin, &source_pos.data()[0], &scal);
}
double fp_emdisk(const std::array<double, 4> f1234, const std::array<double, 2> Photon_motion_constants, const std::array<double, 3> source_pos, const double a_spin, const double mu,
                 const std::array<double, 2> Disk_boundaries){
  double scal = 1.0;
  return __blcoordinate_MOD_pemdisk(f1234.data(), &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &source_pos.data()[1], &source_pos.data()[2], &a_spin, &source_pos.data()[0],
                                    &scal, &mu, &Disk_boundaries.data()[0], &Disk_boundaries.data()[1]);
}

std::array<double, 5> YNOGK(double p, std::array<double, 4> f1234, std::array<double, 2> Photon_motion_constants,std::array<double, 3> source_pos, double a_spin){
  double scal = 1.0;
  double radi, mu, phi, time, sigma;
  __blcoordinate_MOD_ynogk(&p, f1234.data(), &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &source_pos.data()[1],
                                   &source_pos.data()[2], &a_spin, &source_pos.data()[0], &scal, &radi, &mu, &phi, &time, &sigma);
  return {radi, mu, phi, time, sigma};
}

std::array<double, 5> fmetricg(std::array<double, 3> source_position, double a_spin){
  std::array<double, 5> metricg; // expnu, exppsi, expmu1, expmu2, somega
  __blcoordinate_MOD_metricg(&source_position.data()[0], &source_position.data()[1], &source_position.data()[2], &a_spin, &metricg.data()[4],
                             &metricg.data()[0], &metricg.data()[1], &metricg.data()[2], &metricg.data()[3]);
  return metricg;
}

double mu2p(std::array<double, 4> f1234, std::array<double, 2> Photon_motion_constants, std::array<double,3> source_position, double a_spin, int t1, int t2, double mu){
  double scal = 1.0;
  double res = __blcoordinate_MOD_mu2p(&f1234.data()[2], &f1234.data()[1], &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &mu, 
                                       &source_position.data()[1], &source_position.data()[2], &a_spin, &t1, &t2, &scal);
  return res;
}

