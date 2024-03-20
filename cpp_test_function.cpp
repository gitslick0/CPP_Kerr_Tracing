#include <iostream>
#include <cmath>

extern "C" {
  // functions I want to import
  double __blcoordinate_MOD_rms(double *a_spin);
  void __blcoordinate_MOD_initialdirection(double* pr, double* ptheta, double* pphi, double* sinobs, double* muobs, double* a_spin, double* r_obs,
                                           double* velocity, double* lambda, double* q, double* f1234);
}

int main() {
  double rms, a_spin, pr, ptheta, pphi, sinobs, muobs, r_obs, lambda, q;
  double* rms_p = &rms, *a_spin_p = &a_spin, *pr_p = &pr, *ptheta_p = &ptheta, *pphi_p = &pphi, *sinobs_p = &sinobs,
          *muobs_p = &muobs, *r_obs_p = &r_obs, *lambda_p = &lambda, *q_p = &q;
  double velocity[3]; double f1234[4];
  std::cout << "Enter spin factor " << std::endl;
  std::cin >> a_spin;
  rms = __blcoordinate_MOD_rms(a_spin_p);
  // Test if the radius of marginal stability is calculated correctly
  std::cout << "The radius of marginal stability for a= " << a_spin << " is " << rms << std::endl; 

  std::cout << "Enter pr, ptheta, pphi, sinobs, muobs, r_obs, velocity123" << std::endl;
  std::cin >> pr >> ptheta >> pphi >> sinobs >> muobs >> r_obs >> velocity[0] >> velocity[1] >> velocity[2];
  std::cout << "velocity chosen: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;

  __blcoordinate_MOD_initialdirection(pr_p, ptheta_p, pphi_p, sinobs_p, muobs_p, a_spin_p, r_obs_p, velocity, lambda_p, q_p, f1234);

  std::cout << "Initial direction: " << f1234[0] << " " << f1234[1] << " " << f1234[2] << " " << f1234[3] << std::endl;

  return 0;
}
