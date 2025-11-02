/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
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
    zReal = (float *)(aligned_alloc(64, width * sizeof(float)));
    zImag = (float *)(aligned_alloc(64, width * sizeof(float)));
}

LineMandelCalculator::~LineMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    free(zReal);
    free(zImag);
    zReal = NULL;
    zImag = NULL;
    data = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
    // @TODO implement the calculator & return array of integers

    int *pdata = data;
    float *zr = zReal;
    float *zi = zImag;

    for (int i = 0; i < height; i++)
    {
        const float y = y_start + i * dy;

        #pragma omp simd aligned(zr, zi, pdata:64)
        for (int j = 0; j < width; j++)
        {
            zr[j] = 0.0f;
            zi[j] = 0.0f;
            data[i * width + j] = limit; // limit = „zatím neunikl“
        }

        // Vektorizovaná vnitřní smyčka přes sloupce (SIMD-friendly)
        for (int k = 0; k < limit; k++)
        {
            #pragma omp simd aligned(zr, zi, data:64)
            for (int j = 0; j < width; j++)
            {
                // Pokud pixel již divergoval, přeskočíme aktualizaci
                if (pdata[i * width + j] < limit)
                    continue;

                float r2 = zr[j] * zr[j];
                float i2 = zi[j] * zi[j];

                // Pokud právě diverguje, uložíme počet iterací
                if (r2 + i2 > 4.0f)
                {
                    pdata[i * width + j] = k;
                }
                else
                {
                    // Jinak pokračujeme v iteraci
                    float zrTemp = r2 - i2 + (x_start + j * dx);
                    zi[j] = 2.0f * zr[j] * zi[j] + y;
                    zr[j] = zrTemp;
                }
            }
        }
    }

    return data;
}