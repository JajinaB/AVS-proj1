/**
 * @file BatchMandelCalculator.cc
 * @author Jana Brandejsova <xbrand12@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 6. 11. 2025
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"


BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
    BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
    // @TODO allocate & prefill memory
    data = (int *)(aligned_alloc(64, height * width * sizeof(int)));
    zReal = (float *)(aligned_alloc(64, width * sizeof(float)));
    zImag = (float *)(aligned_alloc(64, width * sizeof(float)));
    active = (int *)(aligned_alloc(64, width * sizeof(int)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    data = NULL;
    free(zReal);
    zReal = NULL;
    free(zImag);
    zImag = NULL;
    free(active);
    active = NULL;    
}


int * BatchMandelCalculator::calculateMandelbrot () {
    // @TODO implement the calculator & return array of integers
    int *dataOut = data;
    const int halfHeight = height / 2;

    int tile = 64;

    for (int i = 0; i < halfHeight; ++i) {
        const float y = y_start + i * dy;
        int* row = dataOut + (i * width);
        int* mirrorRow = dataOut + (height - i - 1) * width;

        for (int j = 0; j < width / tile; ++j) {
            int* tileRow = row + (j * tile);

            #pragma omp simd aligned(zReal, zImag, active, tileRow:64)
            for (int x = 0; x < tile; ++x) {
                zReal[x] = 0.0f;
                zImag[x] = 0.0f;
                active[x] = 1;
                tileRow[x] = 0;
            }

            for (int iter = 0; iter < limit; ++iter) {
                int activeCount = 0;

                #pragma omp simd aligned(tileRow, zReal, zImag, active:64) reduction(+:activeCount)
                for (int k = 0; k < tile; ++k) {
                    const float x = x_start + (k + tile * j) * dx;

                    const float r2 = zReal[k] * zReal[k];
                    const float i2 = zImag[k] * zImag[k];

                    const int stillActive = (r2 + i2 <= 4.0f) & active[k];

                    const float zRealTemp = r2 - i2 + x;
                    zImag[k] = 2.0f * zReal[k] * zImag[k] + y;
                    zReal[k] = zRealTemp;

                    tileRow[k] += stillActive;
                    active[k] = stillActive;
                    activeCount += stillActive;
                }

                if (activeCount == 0)
                    break;
            }
        }
        
        #pragma omp simd aligned(row, mirrorRow:64)
        for (int j = 0; j < width; ++j)
            mirrorRow[j] = row[j];
    }

    return pdata;
}