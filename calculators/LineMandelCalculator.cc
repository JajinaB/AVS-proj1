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
}

LineMandelCalculator::~LineMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    data = NULL;    
    free(xVals);
    xVals = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
    // @TODO implement the calculator & return array of integers

    int *pdata = data;
    const halfHeight = height / 2;

    for (int i = 0; i < halfHeight; ++i) {
        const float y = y_start + i * dy;
        int* row = pdata + i * width;
        int* mirrorRow = pdata + (height - i - 1) * width;

        #pragma omp simd aligned(row:64)
        for (int j = 0; j < width; ++j)
            row[j] = limit;

        #pragma omp simd aligned(xVals, row:64)
        for (int j = 0; j < width; ++j) {
            const float x = xVals[j];
            float zr = 0.0f;
            float zi = 0.0f;
            int iter = 0;

            while (iter < limit) {
                const float r2 = zr * zr;
                const float i2 = zi * zi;
                if (r2 + i2 > 4.0f)
                    break;

                const float zrTemp = r2 - i2 + x;
                zi = 2.0f * zr * zi + y;
                zr = zrTemp;
                ++iter;
            }

            row[j] = iter;
            mirrorRow[j] = iter;
        }
    }

    return pdata;
}