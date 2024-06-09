#include <gtest/gtest.h>
#include"./Kerr_Sim/helpers/Classes_lib.h"

TEST(RMSTest, NullTest){
    EXPECT_EQ(frms(0.0), 6.0) << "Calculating minimal stable circular orbit for spin a=0.0, should be 6.0";
    EXPECT_EQ(round(frms(-0.5)*10000)/10000, 7.5546) << "Calculating minimal stable circular orbit for spin a=-0.5, rounded to 4 decimals. Should be 7.5546";
    EXPECT_EQ(round(frms(0.98)*10000)/10000, 1.614) << "Calculating minimal stable circular orbit for spin a=0.98, should be 1.614";
}

TEST(InitialDirectionTest, UpTest){
    EXPECT_EQ(5,5);
}