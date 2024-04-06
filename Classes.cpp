#include"./helpers/Renderer2.h"
#include"./helpers/Classes_lib.h"



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
  arr3d bl_source_position = cart_to_bl({4.0, 0.0, 2.0}, BH.get_a_spin()); // x,y,z coordinates
  std::cout << "Source position " << bl_source_position[0] << ", " << bl_source_position[1] << ", " << bl_source_position[2] << std::endl;
  arr3d source_position = {bl_source_position.data()[0], std::sin(bl_source_position.data()[1]), std::cos(bl_source_position.data()[1])};
  arr3d source_velocity = {0.0, 0.0, 0.0}; // x,y,z velocities
      //
  Emission_Setup ESt = Emission_Setup(source_position, source_velocity, 10, BH.get_a_spin());

  // Check Correctness
  std::cout << "We have " << ESt.get_EmittedPhotons().size() << "Photons" << std::endl;
  ESt.calculate_Photon_rays(1000);
  std::cout << "The first Photon has initial position of " << std::endl;
  ESt.get_EmittedPhotons()[0].displayProperties();
  std::cout << ESt.get_EmittedPhotons()[0].get_total_ray()[1][0] << ", " << ESt.get_EmittedPhotons()[0].get_total_ray()[1][1] << ", " << ESt.get_EmittedPhotons()[0].get_total_ray()[1][2] << std::endl; 

glLineWidth(3.0f);

   const int WIDTH = 800;
   const int HEIGHT = 600;
   const char* TITLE = "OpenGL Window";

   Renderer renderer(WIDTH, HEIGHT, TITLE);
   
   DrawableSphere Sphere = DrawableSphere(BH.get_rHp());
   Sphere.setColor(glm::vec3(0.0, 0.0, 0.0));
   DrawableDisk Disk = DrawableDisk(BH.get_rms(), 100.0, 0.0);
   std::vector<DrawableLine> PhotonDraws = {};
   for(int j = 0; j < ESt.get_EmittedPhotons().size(); j++)
    {
        PhotonDraws.push_back(ESt.get_EmittedPhotons()[j].get_Drawable());
    }
    for(int j = 0; j < PhotonDraws.size(); j++){
        renderer.addDrawable(&PhotonDraws[j]);
    }

   renderer.addDrawable(&Sphere);
   renderer.addDrawable(&Disk);
   //renderer.addDrawable(&Line);

   renderer.render();

  return 0;
}

