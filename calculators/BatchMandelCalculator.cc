/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"


#ifndef BLOCK_H
#define BLOCK_H 32
#endif

#ifndef BLOCK_W
#define BLOCK_W 64
#endif



BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
    BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
    // @TODO allocate & prefill memory
    data = (int *)(aligned_alloc(64, height * width * sizeof(int)));

}

BatchMandelCalculator::~BatchMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    data = NULL;    
}


int * BatchMandelCalculator::calculateMandelbrot () {
    // @TODO implement the calculator & return array of integers
    const float xs = x_start;
    const float ys = y_start;
    const float dx_ = dx;
    const float dy_ = dy;

    for (unsigned bi = 0; bi < height; bi += BLOCK_H)
    {
        for (unsigned bj = 0; bj < width; bj += BLOCK_W)
        {
            for (unsigned i = bi; i < bi + BLOCK_H; ++i)
            {
                float y = ys + i * dy_;
                int *row_out = data + i * width;

                // výpočet počáteční hodnoty x pro daný blok
                float x = xs + bj * dx_;

                // SIMD loop – nyní bez vektorů v paměti, jen aritmetika
                #pragma omp simd aligned(row_out:64)
                for (unsigned j = bj; j < bj + BLOCK_W; ++j)
                {
                    float zr = x;
                    float zi = y;

                    int iter = 0;
                    while (iter < (int)limit)
                    {
                        float r2 = zr * zr;
                        float i2 = zi * zi;
                        if (r2 + i2 > 4.0f) break;

                        float newZi = 2.0f * zr * zi + y;
                        zr = r2 - i2 + x;
                        zi = newZi;

                        iter++;
                    }

                    row_out[j] = iter;
                    x += dx_;        // inkrement místo paměťového přístupu
                }
            }
        }
    }

    return data;
}