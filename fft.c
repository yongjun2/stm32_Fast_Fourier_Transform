/*
 * fft.c
 *
 *  Created on: Jun 9, 2023
 *      Author: jelee
 */



#include <main.h>

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> //realloc
#include <string.h> // memcpy
#include "fft.h"

#define PI 3.1415926535897931160

typedef struct {
  int32_t real;
  int32_t imag;
} Complex;

#define FIXED_POINT_SHIFT 30
#define FIXED_POINT_MULTIPLIER (1 << FIXED_POINT_SHIFT)

int32_t float_to_fixed(float num) {
  return (int32_t)(num * FIXED_POINT_MULTIPLIER);
}

float fixed_to_float(int32_t num) {
  return (float)num / FIXED_POINT_MULTIPLIER;
}


Complex fft_data[IQADC_SAMPLE_SIZE];
float post_fft_data[IQADC_SAMPLE_SIZE];
//float signal_power[41];

int sl, sh;
int thd_idx[8] = {2, 3, 4, 5, 6, 7, 8, 9};
float ssfdr = 0, pwsignal = 0, pwdnoise = 0, pwnoise = 0, pwthd = 0, pwsfdr = 0,
      enob = 0;
int maxidx = 0, maxsfdridx = 0;

#define ARRAYSIZE(A) sizeof(A) / sizeof((A)[0])
float floatArray_max(Complex a[], int size) {
  float max = a[0].real;
  for (int i = 1; i < size; i++)
    if (a[i].real > max)
      max = a[i].real;
  return max;
}


void dit_fft(Complex *buf, int n) {
  int i, j, k;
  int step, block;
  Complex t, u;

  // Bit-reversal permutation
  for (i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j >= bit; bit >>= 1)
      j -= bit;
    j += bit;
    if (i < j) {
      t = buf[i];
      buf[i] = buf[j];
      buf[j] = t;
    }
  }

  // Cooley-Tukey butterfly operations
  for (step = 2; step <= n; step <<= 1) {
    block = step >> 1;
    int32_t angle_fixed = (-2.0 * PI / step) * FIXED_POINT_MULTIPLIER;
    int32_t cosval_fixed = float_to_fixed(cos(fixed_to_float(angle_fixed)));
    int32_t sinval_fixed = float_to_fixed(sin(fixed_to_float(angle_fixed)));

    for (i = 0; i < n; i += step) {
      int32_t w_real = FIXED_POINT_MULTIPLIER;
      int32_t w_imag = 0;

      for (j = 0; j < block; j++) {
        k = i + j;
        t.real = ((int64_t)w_real * buf[k + block].real - (int64_t)w_imag * buf[k + block].imag) >> FIXED_POINT_SHIFT;
        t.imag = ((int64_t)w_real * buf[k + block].imag + (int64_t)w_imag * buf[k + block].real) >> FIXED_POINT_SHIFT;
        u.real = buf[k].real;
        u.imag = buf[k].imag;

        buf[k].real = u.real + t.real;
        buf[k].imag = u.imag + t.imag;
        buf[k + block].real = u.real - t.real;
        buf[k + block].imag = u.imag - t.imag;

        int32_t w_temp_real = ((int64_t)w_real * cosval_fixed - (int64_t)w_imag * sinval_fixed) >> FIXED_POINT_SHIFT;
        w_imag = ((int64_t)w_real * sinval_fixed + (int64_t)w_imag * cosval_fixed) >> FIXED_POINT_SHIFT;
        w_real = w_temp_real;
      }
    }
  }
}

void fft(Complex *buf, int n) {

  // dit_fft(buf, start, end);
  dit_fft(buf, n);
}

void check_time_fft(void)
{
	float sum = 0, average = 0, val = 0, temp = 0;

	for (int i = 0; i < IQADC_SAMPLE_SIZE; i++) {
	  fft_data[i].real = float_to_fixed((float)dump_iq_adc[i]/16384);
	}

	for (int i = 0; i < IQADC_SAMPLE_SIZE; i++) {
		sum = sum + fixed_to_float(fft_data[i].real);
	}
	average = (sum / IQADC_SAMPLE_SIZE);
	//printf("average:%f sum:%f\r\n", average, sum); // ok..

//	for (int i = 0; i < IQADC_SAMPLE_SIZE; i++) {
//		fft_data[i].real = 0;fft_data[i].imag = 0;
//	}

	for (int i = 0; i < IQADC_SAMPLE_SIZE; i++) {
		val = fixed_to_float(fft_data[i].real);
		temp = val - average;
		// fft_data[i].real = (temp * hanning_data[i]);
		fft_data[i].real = float_to_fixed(temp * calculate_hanning[i]);
		fft_data[i].imag = 0.0;
		//printf("idx:%d val:%f temp:%f real:%d\r\n", i, val,temp,  fft_data[i].real);
	  }

//	 for(int i=0; i<20; i++)
//	  printf("index:%d real:%d\r\n", i,fft_data[i].real);


	fft(fft_data, IQADC_SAMPLE_SIZE);

//	printf("after fft process\r\n");
//	for (int i = 20000; i < 20025; i++)
//	    printf("index:%d real:%.7f imag:%.7f\r\n", i,   fixed_to_float(fft_data[i].real),fixed_to_float(fft_data[i].imag));


//  float inputFloat[5] = {3.07434e-09, -2.42417e-08, -9.73708e-08, 9.82207e-08};
//   int32_t fixedValue[5];
//   float floatValue[5];
//
//  for (int i = 0; i < 5; i++) {
//     fixedValue[i] = float_to_fixed(inputFloat[i]);
//     floatValue[i] = fixed_to_float(fixedValue[i]);
//     printf("Fixed Point Value[%d]: %d\n", i, fixedValue[i]);
//     printf("Float Value[%d]: %g\n", i, floatValue[i]);
//   }
//

}
