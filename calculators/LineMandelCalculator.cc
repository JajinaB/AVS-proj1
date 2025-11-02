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
}

LineMandelCalculator::~LineMandelCalculator() {
    // @TODO cleanup the memory
    free(data);
    data = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {
    // @TODO implement the calculator & return array of integers

    int *pdata = data;

    for (int i = 0; i < height; ++i)
    {
        const float y = y_start + i * dy;

        // Vektorizovaná vnitřní smyčka přes sloupce (SIMD-friendly)
        #pragma omp simd aligned(pdata:64)
        for (int j = 0; j < width; ++j)
        {
            const float x = x_start + j * dx;

            // Výpočet Mandelbrotu přímo zde, bez volání mandelbrot()
            float zReal = 0.0f;
            float zImag = 0.0f;
            int iter = 0;

            for (; iter < limit; ++iter)
            {
                const float r2 = zReal * zReal;
                const float i2 = zImag * zImag;

                if (r2 + i2 > 4.0f)
                    break;

                const float zrTemp = r2 - i2 + x;
                zImag = 2.0f * zReal * zImag + y;
                zReal = zrTemp;
            }

            pdata[i * width + j] = iter;
        }
    }

    return data;
}
