/**
 * @file LineMandelCalculator.cc
 * @author Jana Brandejsova <xbrand12@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 3. 11. 2025
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
    BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
    // @TODO allocate & prefill memory
    data = (int *)(aligned_alloc(64, height * width * sizeof(int)));
    xVals = (float *)(aligned_alloc(64, width * sizeof(float)));
    zReal = (float *)(aligned_alloc(64, width * sizeof(float)));
    zImag = (float *)(aligned_alloc(64, width * sizeof(float)));
    active = (int *)(aligned_alloc(64, width * sizeof(int)));

}

LineMandelCalculator::~LineMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    data = NULL;    
    free(xVals);
    xVals = NULL;
    free(zReal);
    zReal = NULL;
    free(zImag);
    zImag = NULL;
    free(active);
    active = NULL;
}

int *LineMandelCalculator::calculateMandelbrot() {
    int *pdata = data;
    const int halfHeight = height / 2;

    for (int i = 0; i < halfHeight; ++i) {
        const float y = y_start + i * dy;
        int* row = pdata + i * width;
        int* mirrorRow = pdata + (height - i - 1) * width;

        #pragma omp simd aligned(zReal, zImag, active, row:64)
        for (int j = 0; j < width; ++j) {
            zReal[j] = 0.0f;
            zImag[j] = 0.0f;
            active[j] = 1;
            row[j] = 0;
        }

        for (int iter = 0; iter < limit; ++iter) {
            int activeCount = 0;

            #pragma omp simd aligned(xVals, row, zReal, zImag, active:64) reduction(+:activeCount)
            for (int j = 0; j < width; ++j) {
                const float x = xVals[j];

                const float r2 = zReal[j] * zReal[j];
                const float i2 = zImag[j] * zImag[j];

                const int stillActive = (r2 + i2 <= 4.0f) & active[j];

                const float zRealTemp = r2 - i2 + x;
                zImag[j] = 2.0f * zReal[j] * zImag[j] + y;
                zReal[j] = zRealTemp;

                row[j] += stillActive;
                active[j] = stillActive;
                activeCount += stillActive;
            }

            if (activeCount == 0)
                break;
        }

        #pragma omp simd aligned(row, mirrorRow:64)
        for (int j = 0; j < width; ++j)
            mirrorRow[j] = row[j];
    }

    return pdata;
}
