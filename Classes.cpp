#include"Classes_lib.h"


int main(int argc, char** argv) {
  //typedefs
  typedef std::array<double, 3> arr3d;
  typedef std::array<double, 4> arr4d;
  typedef std::vector<Photon> vPhoton;
  // initiate random seed
  std::srand(time(NULL));
  // Initialize Emission Setup
      // Initiate BlackHole
  BlackHole BH = BlackHole(0.98);
      // Set up position and Velocity
  arr3d bl_source_position = cart_to_bl({4.0, 0.0, 0.5}, BH.get_a_spin()); // x,y,z coordinates
  arr3d source_position = {bl_source_position.data()[0], std::sin(bl_source_position.data()[1]), std::cos(bl_source_position.data()[1])};
  arr3d source_velocity = {0.0, 0.0, 0.0}; // x,y,z velocities
      //
  Emission_Setup ESt = Emission_Setup(source_position, source_velocity, 1, BH.get_a_spin());

  // Check Correctness
  std::cout << "We have " << ESt.get_EmittedPhotons().size() << "Photons" << std::endl;
  return 0;
}

