#ifndef _Classes_lib_
#define _Classes_lib_

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
  double __pemfinding_MOD_pemfind(const double* f1234, const double* lambda, const double* q, const double* sinobs, const double* muobs, const double* a_spin, const double* robs,
                                 const double* scal, const double* r_in, const double* r_out, const double* muup, const double* mudown, const double* phy1, const double* phy2, const int* caserange,
                                 std::function<double(const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const double* )>* Fp,
                                 const double* paras, const bool* bisection, const double* NN);
}

// Convert fortran wrappers to C++ functions
double frms(const double a_spin){
  return __blcoordinate_MOD_rms(&a_spin);
}
std::array<double, 4> finitial_direction (std::array<double, 3> momenta, std::array<double, 3> source_pos, const double a_spin, std::array<double, 3> source_vel){
  double lambda = 0, q=0;
  std::array<double,4> f1234;
  __blcoordinate_MOD_initialdirection(&momenta.data()[1], &momenta.data()[2], &momenta.data()[0], &source_pos.data()[1], &source_pos.data()[2], &a_spin, &source_pos.data()[0],
                                      source_vel.data(), &lambda, &q, f1234.data());
  return f1234;
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
  double radi = 0.0, mu = 0.0, phi = 0.0, time = 0.0, sigma = 0.0;
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

double pemfind(std::array<double, 4> f1234, std::array<double, 2> Photon_motion_constants, std::array<double,3> source_position, double a_spin, std::array<double, 2> disk_radii, std::array<double, 2> mu_boundaries,
               std::array<double, 2> phi_boundaries, int caserange,
               /*std::function<double(double, std::array<double, 4>, std::array<double, 2>, std::array<double, 3>, double, int, int, std::array<double, 10>)> Fp*/
               std::function<double(const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const double* )> Fp,
                std::array<double, 10> paras, bool bisection, double NN){
                double scal = 1.0;
                return __pemfinding_MOD_pemfind(f1234.data(), &Photon_motion_constants.data()[0], &Photon_motion_constants.data()[1], &source_position.data()[1], &source_position.data()[2],
                &a_spin, &source_position.data()[0], &scal, &disk_radii.data()[0], &disk_radii.data()[1], &mu_boundaries.data()[0], &mu_boundaries.data()[1], &phi_boundaries.data()[0], &phi_boundaries.data()[1],
                &caserange, &Fp, paras.data(), &bisection, &NN);
               }


// Coordinate Related Stuff (mostly obsolete, has been replaced by coordinates.h. Bisection Class might be useful as template later on somewhere)


std::array<double, 4> bl_to_cart(std::array<double,4> coords, double a_spin){
  double x = std::sqrt(coords.data()[0]*coords.data()[0] + a_spin*a_spin) * std::sin(coords.data()[1])*std::cos(coords.data()[2]);
  double y = std::sqrt(coords.data()[0]*coords.data()[0] + a_spin*a_spin) * std::sin(coords.data()[1])*std::sin(coords.data()[2]);
  double z = coords.data()[0] * std::cos(coords.data()[1]);
  double t = coords.data()[3];
  coords = {x,y,z,t};
  return coords;
};

/*
class BisectionMethod1D{
private:
  std::function<double(double)> ThisFunc;
protected:
  void set_ThisFunc(std::function<double(double)> func){
    ThisFunc = func;
  }
public: 
  double lower_bound, upper_bound;
  double estimate;
  double precision = 1.0e-6;
  
  BisectionMethod1D(){};
  BisectionMethod1D(const std::function<double(double)>& inp_func) : ThisFunc(inp_func){
    lower_bound = -10000; upper_bound = 100000; estimate = (lower_bound+upper_bound)/2.0;
  };
  BisectionMethod1D(const std::function<double(double)>& inp_func, double a, double b){
    ThisFunc = inp_func; lower_bound=a; upper_bound=b; estimate = (lower_bound+upper_bound)/2.0;};
  BisectionMethod1D(const std::function<double(double)>& inp_func, std::array<double, 2> bounds){ThisFunc=inp_func; lower_bound=bounds.data()[0]; upper_bound=bounds.data()[1];};
  BisectionMethod1D(const std::function<double(double)>& inp_func, double EP){ThisFunc=inp_func; lower_bound=-10000; upper_bound = 10000; precision = EP;};
  BisectionMethod1D(const std::function<double(double)>& inp_func, double a, double b, double EP){ThisFunc=inp_func; lower_bound=a; upper_bound=b; precision=EP;};

  double bisection(){
    //std::cout << "upper bound at " << upper_bound << " lower bound at " << lower_bound;
    if(ThisFunc(lower_bound)*ThisFunc(upper_bound) >= 0){
      //std::cout << "No sign change detected at lower bound" << lower_bound << " and upper bound " << upper_bound << std::endl;
      return (lower_bound+upper_bound)/2.0;
    }
    double c=lower_bound;
    while(upper_bound-lower_bound >= precision){
      c = (upper_bound+lower_bound)/2.0;
      //std::cout << "Setting center at " << c << std::endl;
      if (std::abs(ThisFunc(c)) <= precision){
        //std::cout << "Function at " << c << " returns " << std::abs(ThisFunc(c)) << std::endl;
        return c;
      }
      else if (ThisFunc(c)*ThisFunc(lower_bound) < 0){
        upper_bound=c;
        //std::cout << "setting upper bound to " << c << std::endl; 
      }
      else {
        //std::cout << "setting lower bound to " << c << std::endl;
        lower_bound = c;
      }
    }
    return c;
  };
};


class BisectionMethod_for_cart_to_bl{
private:
protected:
public: 
  double lower_bound = -0.1, upper_bound = 1.5615;
  double estimate;
  double precision = 1.0e-6;
  double a_spin = 0.0;
  std::array<double, 3> coords = {0.0, 0.0, 0.0};

  BisectionMethod_for_cart_to_bl(){
    estimate = 0.0, lower_bound = -1000.0, upper_bound = 1000.0;
  }
  BisectionMethod_for_cart_to_bl(double a, double b){
    estimate = (a+b)/2.0, upper_bound = a, lower_bound = b;
  }
  BisectionMethod_for_cart_to_bl(double a, double b, double inp_a_spin, std::array<double, 3> inp_coords){
    estimate = (a+b)/2; lower_bound = b, upper_bound = a, coords = coords; a_spin = inp_a_spin;
  }
  BisectionMethod_for_cart_to_bl(std::array<double,3> inp_coords, double inp_a_spin){
    estimate = (lower_bound + upper_bound)/2.0, a_spin = inp_a_spin, coords = inp_coords;
  }
  double root_func(double theta){
    double res = coords[2]*coords[2]*std::tan(theta)*std::tan(theta) + a_spin*a_spin*std::sin(theta)*std::sin(theta) - coords[0]*coords[0] - coords[1]*coords[1];
    return res;
  }

  double bisection(){
    // std::cout << "upper bound at " << upper_bound << " lower bound at " << lower_bound;
    if(root_func(lower_bound)*root_func(upper_bound) >= 0){
      // std::cout << "No sign change detected at lower bound" << lower_bound << " and upper bound " << upper_bound << std::endl;
      return (lower_bound+upper_bound)/2.0;
    }
    double c=lower_bound;
    while(upper_bound-lower_bound >= precision){
      c = (upper_bound+lower_bound)/2.0;
      // std::cout << "Setting center at " << c << std::endl;
      // std::cout << "Function at c is " << root_func(c) << std::endl;
      if (std::abs(root_func(c)) <= precision){
        // std::cout << "Function at " << c << " returns " << std::abs(root_func(c)) << std::endl;
        return c;
      }
      else if (root_func(c)*root_func(lower_bound) < 0){
        upper_bound=c;
        // std::cout << "setting upper bound to " << c << std::endl; 
      }
      else {
        // std::cout << "setting lower bound to " << c << std::endl;
        lower_bound = c;
      }
    }
    return c;
  };

  std::array<double,3> bl_coord(){
    double theta = bisection();
    double radius = coords[2]/std::cos(theta);
    double phi = std::atan(coords[1]/coords[0]);
    return {radius, theta, phi};
  }
};

std::array<double, 3> bl_to_cart(std::array<double, 3> coords, double a_spin){
  double x = std::sqrt(coords[0]*coords[0] + a_spin * a_spin) * std::sin(coords[1])*std::cos(coords[2]);
  double y = std::sqrt(coords[0]*coords[0] + a_spin * a_spin) * std::sin(coords[1])*std::sin(coords[2]);
  double z = coords[0] * std::cos(coords[1]);
  return {x,y,z};
}

std::array<double, 3> cart_to_bl(std::array<double, 3> coords, double a_spin){
  double x = coords[0], y = coords[1], z = coords[2];
  double R1 = std::sqrt(pow(z,4) + 2*pow(y,2)*pow(z,2) + 2*pow(x,2)*pow(z,2) + 2*pow(a_spin,2)*pow(z,2) + pow(y,4) + 2*pow(x,2)*pow(y,2) - 2*pow(a_spin,2)*pow(y,2) + 
                        pow(x,4) - 2*pow(a_spin,2)*pow(x,2) + pow(a_spin,4));
  double radius = std::sqrt(R1 + pow(x,2)+pow(y,2)+pow(z,2)-pow(a_spin,2))/std::sqrt(2);
  double theta = std::acos(z/radius);
  double phi = std::atan(y/x);
  return {radius, theta, phi};
}

std::array<double, 3> cart_to_bl_num(std::array<double,3> coords, double a_spin){
  BisectionMethod_for_cart_to_bl Root_Searcher = BisectionMethod_for_cart_to_bl(coords, a_spin);
  std::array<double,3> new_coords = Root_Searcher.bl_coord();
  return new_coords;
}
*/


// Black Hole class definition
// Input: a_spin = 0; Spin parameter of the black hole
// Methods: get_a_spin(), get_rms, get_rHp, get_rHm, get_rphp, get_rphm
// private properties: _a_spin, _rms, _rHp, _rHm, _rphp, _rphm

double calculate_rho(double radius, double theta, double a_spin){
  return pow((radius*radius + a_spin * a_spin),2) - (radius*radius * a_spin*a_spin + 2*radius*a_spin*a_spin - pow(a_spin, 4)) * std::sin(theta);
}

double calculate_omega(double radius, double theta, double a_spin){
  return 2*a_spin*radius/(calculate_rho(radius, theta, a_spin));
}

double calculate_ut(double radius, double a_spin){
  return (radius * std::sqrt(radius) + a_spin)/(std::sqrt(radius) * std::sqrt(radius*radius - 3*radius + 2*a_spin *std::sqrt(radius)));
}

double calculate_omegadisc(double radius, double a_spin){
  return 1/(radius * std::sqrt(radius) + a_spin);
}

//double energyshift(std::array<double 3> velocity, std::array<double, 3> initial_momenta, std::array<double, 2> Photon_motion_constants, )

double calculate_rms(double a_spin) {
  double Z1 = 1+std::cbrt(1-a_spin*a_spin)*(std::cbrt(1+a_spin) + std::cbrt(1-a_spin));
  double Z2 = sqrt(3*a_spin*a_spin + Z1*Z1);
  double rms = 3+Z2 - copysign(1.0, a_spin)*sqrt((3-Z1)*(3+Z1+2*Z2));
  return rms;
};

double calculate_rHp(double a_spin) {
  return 1+sqrt(1-a_spin*a_spin);
};

double calculate_rHm(double a_spin) {
  return 1-sqrt(1-a_spin*a_spin);
};

double calculate_rphp(double a_spin) {
  return 2*(1+std::cos(2.0/3.0 * std::acos(-a_spin)));
};

double calculate_rphm(double a_spin) {
  return 2*(1+std::cos(2.0/3.0 * std::acos(a_spin)));
};


class AccretionDisk{
private:
  std::array<double, 2> Disk_Boundaries = {0.0, 1000.0};
  std::vector<double> Radius_Bins = {};
public:
  AccretionDisk(std::array<double, 2> input_Boundaries){
    Disk_Boundaries = input_Boundaries;
  }
  AccretionDisk(){
    //pass
  }
  void set_Disk_Boundaries(std::array<double, 2> input_boundaries){
    Disk_Boundaries = input_boundaries;
  }
  void set_Disk_Boundaries(double rin, double rout){
    Disk_Boundaries = {rin, rout};
  }
  void set_Radius_Bins(std::vector<double> input_Bins){
    Radius_Bins = input_Bins;
  }
  std::array<double, 2> get_Disk_Boundaries(){
    return Disk_Boundaries;
  }
  std::vector<double> get_Radius_Bins(){
    return Radius_Bins;
  }
};

class BlackHole{
private:
  double _a_spin, _rms, _rHp, _rHm, _rphp, _rphm;
  AccretionDisk AD = AccretionDisk();

public:
  BlackHole() : _a_spin(0.0), _rms(calculate_rms(0.0)), _rHp(calculate_rHp(0.0)), _rHm(calculate_rHm(0.0)), _rphp(calculate_rphp(0.0)), _rphm(calculate_rphm(0.0)) {};


  BlackHole(double spin) : _a_spin(spin) {
    _rms = calculate_rms(spin); _rHp = calculate_rHp(spin); _rHm = calculate_rHm(spin); _rphp = calculate_rphp(spin); _rphm = calculate_rphm(spin);
    AD = AccretionDisk({_rms, static_cast<double>(INFINITY)});
  };
  BlackHole(double spin, std::array<double,2> Boundaries){
  _rms = calculate_rms(spin); _rHp = calculate_rHp(spin); _rHm = calculate_rHm(spin); _rphp = calculate_rphp(spin); _rphm = calculate_rphm(spin);
  AD = AccretionDisk(Boundaries);
  }

  // Get current spin parameter
  double get_a_spin() const{
    return _a_spin;
  };
  double get_rms() const{
    return _rms;
  };
  double get_rHp() const{
    return _rHp;
  };
  double get_rHm() const{
    return _rHm;
  };
  double get_rphp() const{
    return _rphp;
  };
  double get_rphm() const{
    return _rphm;
  };
  AccretionDisk get_AccretionDisk(){
    return AD;
  }

  void set_a_spin(double new_spin) {
    _a_spin = new_spin; _rms = calculate_rms(new_spin); _rHp = calculate_rHp(new_spin); _rHm = calculate_rHm(new_spin); _rphp=calculate_rphp(new_spin); _rphm=calculate_rphm(new_spin);
  }
  void set_AccretionDisk(AccretionDisk input_Disk){
    AD = input_Disk;
  }
  void set_AccretionDisk(std::array<double,2> input_Boundaries){
    AD = AccretionDisk(input_Boundaries);
  }

  double calculate_frms(double a_spin) {
    return frms(a_spin);
  }
};



// Class Photon
class Photon{
private:
  std::array<double, 2> emission_angles; //{Phi, Theta}
  std::array<double, 3> initial_momenta; //{PPHi, PR, PTheta}
  std::array<double, 4> initial_direction = {0.0, 0.0, 0.0, 0.0}; // Initial direction can later be calculated from the function f1234 once you attach a source to the photon
  std::array<double, 2> motion_constants = {0.0, 0.0};
  double p_total = static_cast<double>(INFINITY);
  double p_emdisk = static_cast<double>(INFINITY);
  double final_destination = -0.5;
  std::vector<std::array<double, 4>> total_ray = {};
  std::vector<CoordVec3> coordinate_ray;
  std::vector<float> time_ray;
  std::vector<glm::vec3> ray_colors;
  DrawableLine Ray;

  
public:
  Photon(double PHI, double THETA){
    double PP = std::cos(THETA); double PR = std::sin(THETA)*std::sin(PHI); double PT = std::sin(THETA)*std::cos(PHI);
    emission_angles = {PHI, THETA};
    initial_momenta = {PP, PR, PT};
  }
  Photon(double PP, double PR, double PT){
    double THETA = std::acos(PP);
    double PHI = 0.0;
    if(std::sin(THETA) != 0.0){
      PHI = std::asin(PR/std::sin(THETA));
    };
    emission_angles = {PHI, THETA};
    initial_momenta = {PP, PR, PT};
  }
  Photon(const std::array<double, 3>& values){
    double THETA = std::acos(values[0]);
    double PHI = 0.0;
    if(std::sin(THETA) != 0.0){
      PHI = std::asin(values[1]/std::sin(THETA));
    };
    emission_angles = {PHI, THETA};
    initial_momenta = values;
  }
  Photon(const std::array<double, 2>& values){
    double PP = std::cos(values[1]); double PR = std::sin(values[1])*std::sin(values[0]); double PT = std::sin(values[1])*std::cos(values[0]);
    emission_angles = values;
    initial_momenta = {PP, PR, PT};
  }
  Photon(){
    /*double RAND1 = static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
    double RAND2 = static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
    double PHI = 2*M_PI * RAND1; "Uniform PHI". Also comment in RAND1 again if you want to use this.
    double PHI = M_PI/2 - M_PI/4 + M_PI/2*RAND1;
    double PHI = M_PI/2; // TOWARDS the black hole
    Different possibilities for PHI */
    double PHI = 0.0;
    //double THETA = std::acos(1-2*RAND2); "uniform Theta". Comment in RAND2 again if you want to use this.
    double THETA = std::acos(1-2*0.5);
    double PP = std::cos(THETA); double PR = std::sin(THETA)*std::sin(PHI); double PT = std::sin(THETA)*std::cos(PHI);
    emission_angles = {PHI, THETA};
    initial_momenta = {PP, PR, PT};
  }
  void set_initial_direction_man(const std::array<double, 4>& values){
    initial_direction = values;
  }
  void set_motion_constants_ind_man(const double lambda, const double q){
    motion_constants = {lambda, q};
  }
  void set_motion_constants_man(const std::array<double, 2> &constants){
    motion_constants = constants;
  }
  void set_p_total(const double total_p){
    p_total = total_p;
  }
  void set_p_emdisk(const double input_pemdisk){
    p_emdisk = input_pemdisk;
  }
  void set_final_destination(const double input_destination){
    final_destination = input_destination;
  }
  void set_total_ray(const std::vector<std::array<double,4>> inp_ray){
    total_ray = inp_ray;
  }
  std::array<double,3> get_initial_momenta(){
    return initial_momenta;
  }
  std::array<double, 2> get_emission_angles(){
    return emission_angles;
  }
  std::array<double, 4> get_initial_direction(){
    return initial_direction;
  }
  std::array<double, 2> get_motion_constants(){
    return motion_constants;
  }
  double get_p_total(){
    return p_total;
  }
  double get_p_emdisk(){
    return p_emdisk;
  }
  double get_final_destination(){
    return final_destination;
  }
  std::vector<std::array<double,4>> get_total_ray(){
    return total_ray;
  }

  void set_coordinate_ray(std::vector<CoordVec3> inp_coordinates){
    this->coordinate_ray = inp_coordinates;
  }
  std::vector<CoordVec3> get_CoordinateRay(){
    return this->coordinate_ray;
  }
  void set_time_ray(std::vector<float> inp_time_ray){
    this->time_ray = inp_time_ray;
  }
  std::vector<float> get_time_ray(){
    return this->time_ray;
  }

  DrawableLine get_Drawable(){
    return Ray;
  }
  void set_Drawable(DrawableLine inp_Line){
    this->Ray = inp_Line;
  }

  std::vector<glm::vec3> get_ray_colors(){
    return this->ray_colors;
  }
  void set_ray_colors(std::vector<glm::vec3> inp_ray_colors){
    this->ray_colors = inp_ray_colors;
  }

  /*void set_initial_direction(const std::array<double, 3>& source_position, const double spin, const std::array<double, 3>& source_velocity){ 
    std::array<double, 4>f1234;  double lambda, q;
    __blcoordinate_MOD_initialdirection(initial_momenta[1], initial_momenta[2], initial_momenta[0], source_position[1], source_position[2], spin, source_position[0], source_velocity, lambda, q, f1234);
    initial_direction = f1234;
  }*/


  void displayProperties() const{
    std::cout << "Emission Angles: { ";
    for (const auto &value : emission_angles){
      std::cout << value << " ";
    }
    std::cout << "}" << std::endl;
    std::cout << "Initial Momenta : { ";
    for (const auto &value : initial_momenta){
      std::cout << value << " ";
    }
    std::cout << "}" << std::endl;
    std::cout << "Photons p_emdisk " << p_emdisk << std::endl;
    std::cout << "Photons p_total " << p_total << std::endl;
    std::cout << "Photons lambda = " << motion_constants[0] << "and q = " << motion_constants[1] << std::endl;
    std::cout << "Photons Initial Direction " << initial_direction[0] << " " << initial_direction[1] << " " << initial_direction[2] << " " << initial_direction[3] << std::endl;
  }
};

// Class Emission Source
class Emission_Source{
private:
  std::array<double, 3> position, velocity; //{robs, sinobs, muobs}, {vr, vtheta, vphi}
  std::vector<Photon> Photons;

public:
  Emission_Source(){
    // Weird that we need an empty/zero constructor in order to declare in Emission_Setup that it has a property of type Emission_Source without
    // initializing immediately....
  }

  Emission_Source(const std::array<double, 3>& input_position, const std::array<double, 3>& input_velocity, int input_numPhotons){
  position = input_position;
  std::cout << "position is set to robs = " << position[0] << "sinobs = " << position[1] << "muobs = " << position[2] << std::endl;
  velocity = input_velocity;
  numPhotons = input_numPhotons;
  for(int ii = 0; ii < numPhotons; ++ii){
      Photon RandomPh = Photon();
      Photons.push_back(RandomPh);
    }
  }

  Emission_Source(const std::array<double, 3>& input_position, const std::array<double, 3>& input_velocity, int input_numPhotons, double a_spin){
    position = input_position;
    std::cout << "position is set to robs = " << position[0] << "sinobs = " << position[1] << "muobs = " << position[2] << std::endl;
    velocity = input_velocity;
    numPhotons = input_numPhotons;
    for (int ii = 0; ii<numPhotons; ++ii){
      Photon RandomPh = Photon();
      RandomPh.set_initial_direction_man(finitial_direction(RandomPh.get_initial_momenta(), position, a_spin, velocity));
      RandomPh.set_motion_constants_man(calculate_photon_motion_constants (RandomPh.get_initial_momenta(), position, a_spin, velocity));
      RandomPh.set_p_total(fp_total(RandomPh.get_initial_direction(), RandomPh.get_motion_constants(), position, a_spin));
      RandomPh.set_p_emdisk(fp_emdisk(RandomPh.get_initial_direction(), RandomPh.get_motion_constants(), position, a_spin, 0.0, {1000.0, calculate_rms(a_spin)}));
      std::cout << "pemdisk" << RandomPh.get_p_emdisk() << std::endl;
      Photons.push_back(RandomPh);
    }
  }
  
  void set_Photons(const std::vector<Photon>& inp_Photons){
    Photons = inp_Photons;
  }

  int numPhotons;
  // Calculate the initial directions of the photons
  void displayProperties() const {
    std::cout << "These are the photons : {" << std::endl;
    int i = 0;
    for (const auto &value : Photons){
      ++i;
      std::cout << "Photon " << i << ": " << std::endl;
      value.displayProperties();
    }
    std::cout << "}" << std::endl;
  }
  void calculate_rays(int n_steps, double a_spin){
    std::vector<Photon> new_Photons = {};
    for (int jj = 0; jj < int(Photons.size()); ++jj){
      Photon curr_Photon = Photons[jj];
      std::vector<std::array<double,4>> curr_Photon_ray = {};
      std::vector<CoordVec3> curr_Photon_coords = {};
      std::vector<float> curr_Photon_time = {};
      std::vector<glm::vec3> curr_Photon_colors = {};
      double p_final = curr_Photon.get_p_total();
      if(curr_Photon.get_p_emdisk() != -1 && curr_Photon.get_p_emdisk() != -2){
        // If we only want the part until it hits the disk
        p_final = curr_Photon.get_p_emdisk();
      }
      for(int kk = 0; kk < n_steps; ++kk){
        double p_curr = p_final/static_cast<double>(n_steps) * kk;
        std::array<double, 5> curr_photon_pos = YNOGK(p_curr, curr_Photon.get_initial_direction(), curr_Photon.get_motion_constants(), position, a_spin);
        std::array<double, 4> cart_pos = bl_to_cart({curr_photon_pos[0], std::acos(curr_photon_pos[1]), curr_photon_pos[2], p_curr}, a_spin);
        CoordVec3 currCoordinate = CoordVec3(cart_pos[0], cart_pos[1], cart_pos[2]);
        glm::vec3 curr_color = glm::vec3(1.0f, 1.0f, 1.0f) * (float)(3.0/(currCoordinate.cart_norm()));
        curr_Photon_ray.push_back(cart_pos);
        curr_Photon_coords.push_back(currCoordinate);
        curr_Photon_time.push_back(curr_photon_pos[3]);
        curr_Photon_colors.push_back(curr_color);

      }
      curr_Photon.set_total_ray(curr_Photon_ray);
      curr_Photon.set_coordinate_ray(curr_Photon_coords);
      curr_Photon.set_time_ray(curr_Photon_time);
      curr_Photon.set_ray_colors(curr_Photon_colors);
      curr_Photon.set_Drawable(DrawableLine(curr_Photon_coords, curr_Photon_time, curr_Photon_colors));
      //curr_Photon.set_Drawable(DrawableLine(curr_Photon_coords, curr_Photon_time));
      new_Photons.push_back(curr_Photon);
    }
  this->set_Photons(new_Photons);
  std::cout << "Changed Photon List!" << std::endl;
  }
  
  std::vector<Photon> get_Photons(){
    return Photons;
  }
  std::array<double, 3> get_position(){
    return position;
  }
  std::array<double, 3> get_velocity(){
    return velocity;
  }
};

class Emission_Setup{
public:
  BlackHole BH;
  Emission_Source ES;
  AccretionDisk AD;

  Emission_Setup(std::array<double, 3>& input_position, const std::array<double, 3>& input_velocity, int input_numPhotons, const double a_spin = 0.0){
    BH = BlackHole(a_spin);
    ES = Emission_Source(input_position, input_velocity, input_numPhotons, a_spin);
    AD = AccretionDisk();
  }

  BlackHole get_BlackHole(){
    return BH;
  }

  Emission_Source get_Emission_Source(){
    return ES;
  }

  AccretionDisk get_AccretionDisk(){
    return AD;
  }
  std::vector<Photon> get_EmittedPhotons(){
    return ES.get_Photons();
  }
  void set_BlackHole(const BlackHole &new_BH){
    BH = new_BH;
  }
  void set_Emission_Source(const Emission_Source &new_ES){
    ES = new_ES;
  }
  void set_AccretionDisk(const AccretionDisk &new_AD){
    AD = new_AD;
  }
  void calculate_Photon_rays(int n_steps){
    ES.calculate_rays(n_steps, BH.get_a_spin());
  }

};


// Drawing Stuff
//
//
//
//
std::vector<std::array<double, 4>> points;
void drawPoints() {
    glClear(GL_COLOR_BUFFER_BIT); // Clear the color buffer

    glBegin(GL_POINTS); // Begin drawing points
    for (const auto& point : points) {
        glColor3f(1.0, 1.0, 1.0); // Set color to white
        glVertex3d(point[0], point[1], point[2]); // Draw the first three coordinates of the point
    }
    glEnd(); // End drawing points
    
   // Draw coordinate axes
    glBegin(GL_LINES); // Begin drawing lines
    // X-axis (red)
    glColor3f(1.0, 0.0, 0.0); // Set color to red
    glVertex3d(-10.0, 0.0, 0.0);
    glVertex3d(10.0, 0.0, 0.0);
    // Y-axis (green)
    glColor3f(0.0, 1.0, 0.0); // Set color to green
    glVertex3d(0.0, -10.0, 0.0);
    glVertex3d(0.0, 10.0, 0.0);
    // Z-axis (blue)
    glColor3f(0.0, 0.0, 1.0); // Set color to blue
    glVertex3d(0.0, 0.0, -10.0);
    glVertex3d(0.0, 0.0, 10.0);
    glEnd(); // End drawing lines
    
   // Draw axis labels
    glPushMatrix(); // Save current transformation matrix
    glColor3f(1.0, 1.0, 1.0); // Set color to white
    glRasterPos3d(10.1, 0.0, 0.0); // Set position for X-axis label
    glutBitmapString(GLUT_BITMAP_HELVETICA_10, (const unsigned char*)"X");
    glRasterPos3d(0.0, 10.1, 0.0); // Set position for Y-axis label
    glutBitmapString(GLUT_BITMAP_HELVETICA_10, (const unsigned char*)"Y");
    glRasterPos3d(0.0, 0.0, 10.1); // Set position for Z-axis label
    glutBitmapString(GLUT_BITMAP_HELVETICA_10, (const unsigned char*)"Z");
    glPopMatrix(); // Restore previous transformation matrix 

    glutSwapBuffers(); // Swap the front and back buffers to display the rendered image
}

void init_canvas() {
    glClearColor(0.0, 0.0, 0.0, 0.0); // Set the clear color (black)
    glPointSize(5.0); // Set point size to 5 pixels
    glMatrixMode(GL_PROJECTION); // Set the matrix mode to projection
    glLoadIdentity(); // Load the identity matrix
    gluPerspective(45.0, 1.0, 1.0, 100.0); // Set perspective projection
    glMatrixMode(GL_MODELVIEW); // Set the matrix mode to modelview
    glLoadIdentity(); // Load the identity matrix
    gluLookAt(0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); // Set the camera position and orientation
};



#endif