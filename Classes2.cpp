#include"Classes.h"


int main(int argc, char** argv) {
  std::srand(time(NULL));
  double a_spin = 0.98;
  BlackHole BH = BlackHole(a_spin);
  double a = BH.get_a_spin();
  double rHp = BH.get_rHp();
  std::cout << "a = " << BH.get_a_spin() << ", rms = " << BH.get_rms() << ", rHp = " << rHp << ", rHm = " << BH.get_rHm() << std::endl;
  /*std::cout << "The spin parameter is set to " << a << std::endl;
  std::cout << "The outer event horizon is at " << rHp << std::endl;
  std::cout << "The inner event horizon is at " << BH.get_rHm() << std::endl;
  std::cout << "The radius of marginal stability is at " << BH.get_rms() << std::endl;
  std::cout << "The radius of bound photon orbits is at " << BH.get_rphp() << " and " << BH.get_rphm() << std::endl; */
  std::array<double, 3> source_position_cart = {5.0, 0.0, 2.0};
  std::array<double, 3> source_position_bl = cart_to_bl(source_position_cart, a_spin);
  std::array<double, 3> source_pos = {source_position_bl.data()[0], std::sin(source_position_bl.data()[1]), std::cos(source_position_bl.data()[2])};
  std::array<double, 5> metricg = fmetricg(source_pos, a_spin);
  std::cout << "metricg: " << " expnu = " << metricg.data()[0] << " exppsi = " << metricg.data()[1] << " expnu1 = " << metricg.data()[2] <<
    " expnu2 = " << metricg.data()[3] << " somega = " << metricg.data()[4] << std::endl; 

  std::array<double, 3> source_vecocity = {0.0, 0.0, 0.0};
  std::array<double,3> source_position = {source_position_bl.data()[0], std::sin(source_position_bl.data()[1]), std::cos(source_position_bl.data()[1])};
  Emission_Source ES = Emission_Source(source_position, source_vecocity, 1, a_spin);
  ES.displayProperties();
  Photon Da_Photon = ES.get_Photons().data()[0];
  double p_disk = mu2p(Da_Photon.get_initial_direction(), Da_Photon.get_motion_constants(), ES.get_position(), a_spin, 0, 0, 0.3);
  std::cout << "p_disk = " << mu2p << std::endl;

  return 0;
}

