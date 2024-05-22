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
    CoordVec3 bl_source_position_cart = CoordVec3(2.5, 0.0, 0.001);
    CoordVec3 bl_source_position = bl_source_position_cart.to_bl(BH.get_a_spin());
    //CoordVec3 bl_source_position2 = CoordVec3(0.0, 0.0, 4.5).to_bl(BH.get_a_spin());
    //arr3d bl_source_position = cart_to_bl({10.0, 0.0, 4.5}, BH.get_a_spin()); // x,y,z coordinates
    //std::cout << "Source position " << bl_source_position[0] << ", " << bl_source_position[1] << ", " << bl_source_position[2] << std::endl;
    arr3d source_position = {bl_source_position.x, std::sin(bl_source_position.y), std::cos(bl_source_position.y)};
    arr3d source_velocity = {0.0, 0.0, 0.0}; // x,y,z velocities
    //arr3d source_position2 = {bl_source_position2.x, std::sin(bl_source_position2.y), std::cos(bl_source_position2.y)};
      //
    Emission_Setup ESt = Emission_Setup(source_position, source_velocity, 200, BH.get_a_spin());
    //Emission_Setup ESt2 = Emission_Setup(source_position2, source_velocity, 100, BH.get_a_spin());

    // Calculate Photon Rays
    ESt.calculate_Photon_rays(400);
    //ESt2.calculate_Photon_rays(500);

    glLineWidth(3.0f);

    const int WIDTH = 800;
    const int HEIGHT = 600;
    const char* TITLE = "OpenGL Window";

    Renderer renderer(WIDTH, HEIGHT, TITLE);
   
    DrawableSphere Sphere = DrawableSphere(BH.get_rHp());
    Sphere.setColor(glm::vec3(0.0, 0.0, 0.0));
    DrawableDisk Disk = DrawableDisk(BH.get_rms(), 100.0, 0.0);
    DrawableBLSphere Draw_BH = DrawableBLSphere(BH.get_rHp(), BH.get_a_spin());
    std::vector<DrawableLine> PhotonDraws = {};
    //std::vector<DrawableLine> PhotonDraws2 = {};
    for(int j = 0; j < ESt.get_EmittedPhotons().size(); j++)
    {
        PhotonDraws.push_back(ESt.get_EmittedPhotons()[j].get_Drawable());
        //PhotonDraws2.push_back(ESt2.get_EmittedPhotons()[j].get_Drawable());
    }
    for(int j = 0; j < PhotonDraws.size(); j++){
        renderer.addDrawable(&PhotonDraws[j]);
        //renderer.addDrawable(&PhotonDraws2[j]);
    }

    renderer.addDrawable(&Sphere);
    //renderer.addDrawable(&Disk);
    renderer.addDrawable(&Draw_BH);
    //renderer.addDrawable(&Line);

    std::cout << "Starting to render" << std::endl;

    renderer.render();

  return 0;
}

