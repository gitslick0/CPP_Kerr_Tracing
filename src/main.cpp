#include"./Kerr_Sim/helpers/Renderer2.h"
#include"./Kerr_Sim/helpers/Classes_lib.h"



int main() {
    //typedefs
    typedef std::array<double, 3> arr3d;

    int n_sim = 200; //number of steps simulated along each Photons trajectory

    // initiate random seed
    std::srand(time(NULL));
    
    // Initialize Emission Setup
      // Initiate BlackHole
    std::cout << "Set the Spin of the Black Hole between (-0.998, 0.998)" << std::endl;
    float input_a_spin;
    std::cin >> input_a_spin; 
    BlackHole BH = BlackHole(input_a_spin);

      // Set up position and Velocity of source
    std::cout << "Set the (cartesian) position of the source (float x,y,z)\n note that the point must lie outside the event horizon for the visualization to work properly" << std::endl;
    float input_x, input_y, input_z;
    std::cin >> input_x >> input_y >> input_z;
    
    CoordVec3 bl_source_position_cart = CoordVec3(input_x, input_y, input_z);
    CoordVec3 bl_source_position = bl_source_position_cart.to_bl(BH.get_a_spin());

    arr3d source_position = {bl_source_position.x, std::sin(bl_source_position.y), std::cos(bl_source_position.y)};
    arr3d source_velocity = {0.0, 0.0, 0.0}; // x,y,z velocities

    Emission_Setup ESt = Emission_Setup(source_position, source_velocity, 200, BH.get_a_spin());


    std::cout << "int: How many Photons would you like to simulate? (Note: Visualization can only handle ~ 1000 Photons (depending on your specs))" << std::endl;
    int num_Photons;
    std::cin >> num_Photons;

    // Calculate Photon Rays
    ESt.calculate_Photon_rays(num_Photons);
    //ESt2.calculate_Photon_rays(500);

    glLineWidth(3.0f);

    const int WIDTH = 800;
    const int HEIGHT = 600;
    const char* TITLE = "BlackHole ??";

    Renderer renderer(WIDTH, HEIGHT, TITLE);
    DrawableBLSphere Draw_BH = DrawableBLSphere(BH.get_rHp(), BH.get_a_spin());

    
    std::vector<DrawableLine> PhotonDraws = {};
    for(int j = 0; j < int(ESt.get_EmittedPhotons().size()); j++)
    {
        PhotonDraws.push_back(ESt.get_EmittedPhotons()[j].get_Drawable());
        //PhotonDraws2.push_back(ESt2.get_EmittedPhotons()[j].get_Drawable());
    }
    for(int j = 0; j < int(PhotonDraws.size()); j++){
        renderer.addDrawable(&PhotonDraws[j]);
        //renderer.addDrawable(&PhotonDraws2[j]);
    }

    renderer.addDrawable(&Draw_BH);

#ifdef USE_DISK
    DrawableDisk Disk = DrawableDisk(BH.get_rms(), 100.0, 0.0);
    renderer.addDrawable(&Disk);
#endif

    std::cout << "Starting to render" << std::endl;

    renderer.render();

  return 0;
}

