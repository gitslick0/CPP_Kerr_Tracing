# include <iostream>
#include <fitsio.h>
#include <string.h>
#include <cstdio>
#include <gtest/gtest.h>

int create_test_fits(std::string filename){
    fitsfile *fptr;   // Pointer to the FITS file; defined in fitsio.h
    int status = 0;   // CFITSIO status value MUST be initialized to zero!
    long naxes[2] = {300, 200};  // Image dimensions
    long fpixel[2] = {1, 1};     // First pixel to write
    int bitpix = SHORT_IMG;      // 16-bit unsigned short pixels
    int nelements = naxes[0] * naxes[1];  // Number of pixels in the image
    short *array;     // Dynamically allocate memory for the image array

    // Allocate memory for the image
    array = new short[nelements];
    for (int i = 0; i < nelements; i++) {
        array[i] = i % 1000;  // Fill the array with some dummy data
    }

    // Create a new FITS file
    if (fits_create_file(&fptr, filename.c_str(), &status)) {
        fits_report_error(stderr, status); // Print error message
        return status;
    }

    // Create the primary array image (2D short integer image)
    if (fits_create_img(fptr, bitpix, 2, naxes, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Write the array of unsigned short integers to the image
    if (fits_write_pix(fptr, TSHORT, fpixel, nelements, array, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Close the FITS file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
        return status;
    }

    // Free the allocated memory
    delete[] array;
    
    return 0;
};

TEST(CREATE_TESTFILE, Test1){
    // Create and remove a Testfile.fits
    // Create

    int cr_res = create_test_fits("Testfile.fits");
    EXPECT_EQ(cr_res, 0);
    
    // remove
    cr_res = remove("Testfile.fits");
    EXPECT_EQ(cr_res, 0);
}