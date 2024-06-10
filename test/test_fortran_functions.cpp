#include <gtest/gtest.h>
#include"./Kerr_Sim/helpers/Classes_lib.h"
#include <time.h>

TEST(RMSTest, NullTest){
    EXPECT_EQ(frms(0.0), 6.0) << "Calculating minimal stable circular orbit for spin a=0.0, should be 6.0";
    EXPECT_EQ(round(frms(-0.5)*10000)/10000, 7.5546) << "Calculating minimal stable circular orbit for spin a=-0.5, rounded to 4 decimals. Should be 7.5546";
    EXPECT_EQ(round(frms(0.98)*10000)/10000, 1.614) << "Calculating minimal stable circular orbit for spin a=0.98, should be 1.614";
}

TEST(Emission_Setup_and_ZAM_Photons_Test, UpTest){
    // set random seed
    srand(time(NULL));

    // get random spin parameter
    double a_spin = (static_cast<double>(std::rand())/static_cast<double>(RAND_MAX) -0.5) * 0.998;

    // get random source robs:
    double rand_robs = static_cast<double>(std::rand())/static_cast<double>(RAND_MAX) * 15 + 2.5; // robs between 2.5 and 17.5
    double theta = 0.0; // set source on axis of symmetry

    // Initiate Emision_Setup
    Emission_Source ES = Emission_Source({rand_robs, std::sin(theta), std::cos(theta)}, {0.0, 0.0, 0.0}, 300, a_spin);

    // Test if all Photons have zero angular momentum (bc. they are on axis of symmetry)
    for (auto &Ph : ES.get_Photons()){
        EXPECT_EQ(Ph.get_motion_constants().data()[0], 0.0);
    }
}