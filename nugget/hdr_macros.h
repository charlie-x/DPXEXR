
////////////////////////////////////////////////////////////////////////////////////// 
// License Terms for Academy Color Encoding System Components                       //
//                                                                                  //
// Academy Color Encoding System (ACES) software and tools are provided by the      //
// Academy under the following terms and conditions: A worldwide, royalty-free,     //
// non-exclusive right to copy, modify, create derivatives, and use, in source and  //
// binary forms, is hereby granted, subject to acceptance of this license.          //
// Performance of any of the aforementioned acts indicates acceptance to be bound   //
// by the following terms and conditions:                                           //
//                                                                                  //
// Copies of source code, in whole or in part, must retain the above copyright      //
// notice, this list of conditions and the Disclaimer of Warranty.                  //
//                                                                                  //
// Use in binary form must retain the copyright notice (below), this list of        //
// conditions and the Disclaimer of Warranty in the documentation and/or other      //
// materials provided with the distribution.                                        //
//                                                                                  //
// * Nothing in this license shall be deemed to grant any rights to trademarks,     //
// copyrights, patents, trade secrets or any other intellectual property of         //
// A.M.P.A.S. or any contributors, except as expressly stated herein.               //
//                                                                                  //
// * Neither the name "A.M.P.A.S." nor the name of any other contributors to this   //
// software may be used to endorse or promote products derivative of or based on    //
// this software without express prior written permission of A.M.P.A.S. or the      //
// contributors, as appropriate.                                                    //
//                                                                                  //
// * This license shall be construed pursuant to the laws of the State of           //
// California, and any disputes related thereto shall be subject to the             //
// jurisdiction of the courts therein.                                              //
//                                                                                  //
// Copyright © 2012 Academy of Motion Picture Arts and Sciences (A.M.P.A.S.).       //
// Portions contributed by others as indicated. All rights reserved.                //
//                                                                                  //
// Disclaimer of Warranty: THIS SOFTWARE IS PROVIDED BY A.M.P.A.S. AND CONTRIBUTORS //
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,    //
// THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND //
// NON-INFRINGEMENT ARE DISCLAIMED. IN NO EVENT SHALL A.M.P.A.S., OR ANY            //
// CONTRIBUTORS OR DISTRIBUTORS, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    //
// SPECIAL, EXEMPLARY, RESITUTIONARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT  //
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR   //
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF           //
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE  //
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF        //
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                       //
//                                                                                  //
// WITHOUT LIMITING THE GENERALITY OF THE FOREGOING, THE ACADEMY SPECIFICALLY       //
// DISCLAIMS ANY REPRESENTATIONS OR WARRANTIES WHATSOEVER RELATED TO PATENT OR      //
// OTHER INTELLECTUAL PROPERTY RIGHTS IN THE ACADEMY COLOR ENCODING SYSTEM, OR      //
// APPLICATIONS THEREOF, HELD BY PARTIES OTHER THAN A.M.P.A.S.,WHETHER DISCLOSED OR //
// UNDISCLOSED.                                                                     //
//////////////////////////////////////////////////////////////////////////////////////  


/////////////////////////////////////////////////////////////////////////////
//
//
// Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, AMPAS
// Author: Gary Demos 2010,2011,2012,2013,2014, 2015
//
// Algorithms and code herein also are copyright Gary Demos, and some portions may retain a private use and copyright
//
// Further developed in collaboration with the ASC Technology Committee, beginning March 2015
//
// Note: patents may apply to ideas contained herein.
// Such patents may include patents issued and/or pending by Gary Demos.
// No patent licenses are granted nor implied due to open availability of this code.
//
//
//////////////////////////////////////////////////////////////////////////////


/***********************************************************************************************************/
/* select one of the following: */
#define BYPASS_LMT /* this selects solely the rendering nugget */
//#define SIMPLE_LMT
//#define GAMMA_AND_MAT
//#define MODERATE_LMT
//#define FULL_LMT
//#define DOUBLE_LMT


#define BRIGHT_HIGHLIGHTS_IN_TONE_CURVE /* if defined: use bright highlights, if not defined: use lower highlights */
//#define FACES_HIGHLIGHTS /* if defined: desaturate highlights for faces, if not defined: more color saturation in bright colors (does not apply to BYPASS_LMT nor GAMMA_AND_MAT) */

#define RADIOMETRIC_ODTS /* RGB ratio-preserving ODTs (chromaticity preserving) */
//#define PQ_GAMMA /* use PQ gamma HDR instead of gamma exponent */

//***** NOTE: with some GPUs when using an OpenCL version of this code, it is necessary to make the parent .cl file (that includes this file) continuously different in order to force a re-interpret.
// Any previous matching version will fall back to the corresponding cache  Otherwise no change will be observed.
// Recommended is to add a letter or number to the end of a comment line in the ".cl" file each time this file is modified


/***********************************************************************************************************/

#define value_of_pi 3.1415926535f


#define sin_function_180deg(a, b, in_val, slow_val) \
 if (in_val > b) { \
   slow_val = 1.0; \
 } else { \
   if (in_val < a) { \
     slow_val = 0.0; \
   } else { \
     slow_val = 0.5 * (1.0 + sinf(-0.5 * value_of_pi + value_of_pi * (in_val - a)/(b-a))); \
   } \
 }

#define sin_function_360deg(a, b, in_val, slow_val) \
 if (in_val > b) { \
   slow_val = 0.0; \
 } else { \
   if (in_val < a) { \
     slow_val = 0.0; \
   } else { \
     slow_val = 0.5 * (1.0 + sinf(-0.5 * value_of_pi + 2.0 * value_of_pi * (in_val - a)/(b-a))); \
   } \
 }

#define COLOR_MATRIX_TIMES_VEC(rgb_out, MATX, rgb_in) \
  { \
   float MTX_B[9], MTX_C[9], rgb_masked[3]; \
   for (i=0; i<9; i++) { \
/* note: it is expected that diagonal matrix terms will always be positive, but this is implemented as a general operator anyway */ \
     if (MATX[i] >= 0.0) { \
       MTX_B[i] = MATX[i]; \
       MTX_C[i] = 0.0; \
     } else { /* MATX[i] < 0.0 */ \
       MTX_B[i] = 0.0; \
       MTX_C[i] = MATX[i]; /* negative terms go into MTX_C */ \
     } /* MATX[i] >= 0 or not */ \
   } /* i */ \
   for (i=0; i<3; i++) { \
     if (rgb_in[i] >= 0.0) { \
       rgb_masked[i] = rgb_in[i]; \
     } else { /* rgb_in[i] negative */ \
       rgb_masked[i] = 0.0; /* mask out negative numbers for negative matrix terms (they fold over positive) */ \
     } /* rgb_in[i] >= 0.0 or not */ \
   } /* i */ \
   rgb_out[0] = MTX_B[0] * rgb_in[0] +  MTX_B[1] * rgb_in[1] +  MTX_B[2] * rgb_in[2]; \
   rgb_out[1] = MTX_B[3] * rgb_in[0] +  MTX_B[4] * rgb_in[1] +  MTX_B[5] * rgb_in[2]; \
   rgb_out[2] = MTX_B[6] * rgb_in[0] +  MTX_B[7] * rgb_in[1] +  MTX_B[8] * rgb_in[2]; \
/* now add in negative matrix terms times positively-masked (non-negative, with negative rgb values set to zero) rgb values */ \
   rgb_out[0] = rgb_out[0] + MTX_C[0] * rgb_masked[0] + MTX_C[1] * rgb_masked[1] + MTX_C[2] * rgb_masked[2]; \
   rgb_out[1] = rgb_out[1] + MTX_C[3] * rgb_masked[0] + MTX_C[4] * rgb_masked[1] + MTX_C[5] * rgb_masked[2]; \
   rgb_out[2] = rgb_out[2] + MTX_C[6] * rgb_masked[0] + MTX_C[7] * rgb_masked[1] + MTX_C[8] * rgb_masked[2]; \
  }



/*************************************************************************************************************************/
/* compute slow_in_out to eliminate the possibility of visible mach bands (first derivative discontinuities) across neutral (when crossing between hue regions as defined by which r,g, or b primary is largest) */ \
/* goes from zero at -pi/2 for min/max = 1.0 to one at +pi/2 for min/max = (1.0-.2) or less */ \
/* this version uses the second order two-term norm instead of the minimum of colmin1 and colmin2 (also eliminates first derivative discontinuities at the crossover of min of colmin1 vs colmin2) */
#define compute_two_min_norm(colmin1, colmin2, two_min_norm) \
{ \
        float sign_of_colmin1, sign_of_colmin2, sum_of_two_squares, sign_of_sum_of_two_squares; \
        sign_of_colmin1 = 1.0; if (colmin1 < 0.0) { sign_of_colmin1 = -1.0; } \
        sign_of_colmin2 = 1.0; if (colmin2 < 0.0) { sign_of_colmin2 = -1.0; } \
        sum_of_two_squares = sign_of_colmin1 * colmin1 * colmin1 + sign_of_colmin2 * colmin2 * colmin2;\
        sign_of_sum_of_two_squares = 1.0; if (sum_of_two_squares < 0.0) { sign_of_sum_of_two_squares = -1.0; } \
        two_min_norm = sign_of_sum_of_two_squares * powf(sign_of_sum_of_two_squares * sum_of_two_squares * 0.5, 0.5); /* sqrt((a^2+b^2)/2) with negative processing */ \
}

#define compute_slow_in_out_b(sio_depth, colmin1, colmin2, colmax) \
{ \
        float two_min_norm; \
        compute_two_min_norm(colmin1, colmin2, two_min_norm) \
        slow_in_out = 0.5 + 0.5 * sinf(value_of_pi * (1.0 - (MAX(0.0, MIN(1.0, ( two_min_norm/MAX(1e-6, colmax) - (1.0 - sio_depth) ) / sio_depth)))) - value_of_pi/2.0); \
}

#define compute_slow_in_out_c(sio_depth, colmin1, colmin2, colmax) \
{ \
        float two_min_norm; \
        compute_two_min_norm(1.0 / MAX(1e-12, colmin1), 1.0 / MAX(1e-12, colmin2), two_min_norm) /* use inverted minima to make r, g, and b have similar saturation to y, c, and m */ \
        two_min_norm = 1.0 / MAX(1e-12, two_min_norm); /* restore inversion */ \
        slow_in_out = 0.5 + 0.5 * sinf(value_of_pi * (1.0 - (MAX(0.0, MIN(1.0, ( two_min_norm/MAX(1e-6, colmax) - (1.0 - sio_depth) ) / sio_depth)))) - value_of_pi/2.0); \
}

#define compute_slow_in_out_d(sio_depth, colmin, colmax) \
        slow_in_out = 0.5 + 0.5 * sinf(value_of_pi * (1.0 - (MAX(0.0, MIN(1.0, ( colmin/MAX(1e-9, colmax) - (1.0 - sio_depth) ) / sio_depth)))) - value_of_pi/2.0);

/* compute slow_in_out to eliminate the possibility of visible mach bands (first derivative discontinuities) */ \
/* goes from zero at -pi/2 for min/max = 1.0 to one at +3pi/2 for min/max = (1.0-sio_depth) or less */

#define compute_slow_in_out_360_b(sio_depth, colmin1, colmin2, colmax) \
{ \
        float two_min_norm; \
        compute_two_min_norm(colmin1, colmin2, two_min_norm) \
        slow_in_out_360 = 0.5 + 0.5 * sinf(2.0 * value_of_pi * (1.0 - (MAX(0.0, MIN(1.0, ( two_min_norm/MAX(1e-6, colmax) - (1.0 - sio_depth) ) / sio_depth)))) - value_of_pi/2.0); \
}

/***********************************************************************************************************/

#define process_hue_region_matrix(sat_boost_amount, darken_sat, sat_bias, strength) \
              float off_diag_sat = strength * sat_boost_amount + sat_bias; \
              float diag_sat = 1.0 - 2.0 * off_diag_sat; \
/* darken or brighten saturated colors */ \
              off_diag_sat = (1.0 - strength * darken_sat) * off_diag_sat; \
              diag_sat     = (1.0 - strength * darken_sat) * diag_sat; \
\
              float SAT_MTX[9], RGB_ADJUSTED[3]; \
              SAT_MTX[0] = diag_sat;            SAT_MTX[1] = off_diag_sat; SAT_MTX[2] = off_diag_sat; \
              SAT_MTX[3] = off_diag_sat;        SAT_MTX[4] = diag_sat;     SAT_MTX[5] = off_diag_sat; \
              SAT_MTX[6] = off_diag_sat;        SAT_MTX[7] = off_diag_sat; SAT_MTX[8] = diag_sat; \
              COLOR_MATRIX_TIMES_VEC(RGB_ADJUSTED, SAT_MTX, RGB)

#define process_sat_matrix(desat_amount, rgb_desat) \
 { \
              float off_diag_sat = desat_amount; \
              float diag_sat = 1.0 - 2.0 * off_diag_sat; \
\
              float SAT_MTX[9]; \
              SAT_MTX[0] = diag_sat;            SAT_MTX[1] = off_diag_sat; SAT_MTX[2] = off_diag_sat; \
              SAT_MTX[3] = off_diag_sat;        SAT_MTX[4] = diag_sat;     SAT_MTX[5] = off_diag_sat; \
              SAT_MTX[6] = off_diag_sat;        SAT_MTX[7] = off_diag_sat; SAT_MTX[8] = diag_sat; \
              COLOR_MATRIX_TIMES_VEC(rgb_desat, SAT_MTX, RGB) \
 }



/***********************************************************************************************************/
#define MATMUL(rgb_out, MATX, rgb_in) \
   rgb_out[0] = MATX[0] * rgb_in[0] +  MATX[1] * rgb_in[1] +  MATX[2] * rgb_in[2]; \
   rgb_out[1] = MATX[3] * rgb_in[0] +  MATX[4] * rgb_in[1] +  MATX[5] * rgb_in[2]; \
   rgb_out[2] = MATX[6] * rgb_in[0] +  MATX[7] * rgb_in[1] +  MATX[8] * rgb_in[2];


/************************************************************************************************************/
#define TINY 1e-12
#define LARGER_TINY 1e-9

/* HDR Parameters ********************************************************************************************************/
#define GAMMA_BOOST_HDR 1.07 /* compensate for asymptotic function, should not use values below 1. */
#define OFF_DIAG_HDR -.005 /* very sensitive, negative numbers increase color saturation, positive numbers desaturate */
#define HALF_WAY_HDR 2.0 /* point at which asymptotic function for highlights goes to 1/2 */
#define HDR_GAIN .41 /* scale factor, adjust depending on the dynamic range of the display/projector, and the value of the gamma boost above */


/* MDR Parameters ********************************************************************************************************/
#define GAMMA_BOOST_MDR 1.045 /* compensate for asymptotic function, should not use values below 1.0 */
#define OFF_DIAG_MDR -.02 /* very sensitive, negative numbers increase color saturation, positive numbers desaturate */
#define HALF_WAY_MDR 1.4 /* point at which asymptotic function for highlights goes to 1/2 */
#define MDR_GAIN .7 /* scale factor, adjust depending on the dynamic range of the display/projector, and the value of the gamma boost above */

#ifdef PQ_GAMMA
/* PQ 2084 */
/* Base functions from SMPTE ST 2084-2014, aka PQ_gamma Curve */
/* Cribbed from one of the ACES 1.0 HDR ODT's */

#define pq_m1 0.1593017578125 /* ( 2610.0 / 4096.0 ) / 4.0 */
#define pq_m2 78.84375 /* ( 2523.0 / 4096.0 ) * 128.0 */
#define pq_c1 0.8359375 /* 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0 */
#define pq_c2 18.8515625 /* ( 2413.0 / 4096.0 ) * 32.0 */
#define pq_c3 18.6875 /* ( 2392.0 / 4096.0 ) * 32.0 */
#define pq_C 100.0 /* 100.0 / 10000.0 */

/* Converts from linear cd/m^2 (assuming 10 = 0.1 = outputLAD to the non-linear perceptually quantized space */
/* Note that this is in float, and assumes normalization from 0 - 1 */
/* (0 - pq_C for linear) and does not handle the integer coding in the Annex sections of SMPTE ST 2084-2014 */
/* Note that this does NOT handle any of the signal range considerations from SMPTE 2084 - this returns full range (0 - 1) */

 #define APPLY_GAMMA(PQ_Out, PQ_In) \
 { \
   float L = PQ_In / pq_C; \
   float Lm = powf( L, pq_m1 ); \
   float N = ( pq_c1 + pq_c2 * Lm ) / ( 1.0 + pq_c3 * Lm ); \
   PQ_Out = powf( N, pq_m2 ); \
   if (PQ_Out > 1.0) { PQ_Out = 1.0; } \
 }

 #define process_gamma(ignore_the_argument) \
 { \
     APPLY_GAMMA(rOut, RGB_VEC[0]) \
     APPLY_GAMMA(gOut, RGB_VEC[1]) \
     APPLY_GAMMA(bOut, RGB_VEC[2]) \
 }


#else /* exponent gamma */

 #define HDR_DISPLAY_GAMMA 2.4
 #define MDR_DISPLAY_GAMMA 2.4

 #define process_gamma(gamma_exponent_for_display) \
 { \
     rOut = powf(RGB_VEC[0], 1.0/gamma_exponent_for_display); \
     gOut = powf(RGB_VEC[1], 1.0/gamma_exponent_for_display); \
     bOut = powf(RGB_VEC[2], 1.0/gamma_exponent_for_display); \
 }

#endif /* PQ_GAMMA or not */


#ifdef RADIOMETRIC_ODTS

/* a radiometric proportion of 1.0 is fully RGB-relative-ratio preserving, and thus is chromaticity preserving */
/* a radiometric proportion of 0.0 fully desaturates bright colors */
#define MDR_RADIOMETRIC_PROPORTION 0.6 /* darker displays can have more color saturation in bright colors */
#define HDR_RADIOMETRIC_PROPORTION 1.0//0.4 /* brighter displays appear more "colorful" (the Hunt effect), and thus bright colors can be somewhat less saturated */


#define ROOM_BRIGHTENING 0.0 /* a room_brightening of 0.0 is for dark surround viewing, 1.0 is for viewing in a bright room (applied to radiometric portion) */

//#define NORM_CORRECTION  .9 /* found empirically, scale factor corresponding to inverse of 1.111, the practical norm maximum */

/***********************************************************************************************************/
#define xy_preserving_odt(gain, half_way, non_radiometric_asymp_gamma_boost, room_brightening, aces_to_device_mat, \
 display_gamma, radiometric_proportion) \
\
/* input is rgbOut[3], output is rOut, gOut, bOut */ \
\
        { \
\
          if (rgbOut[0]<0.0) {rgbOut[0]=0.0;} /* clip values up to zero */ \
          if (rgbOut[1]<0.0) {rgbOut[1]=0.0;} /* clip values up to zero */ \
          if (rgbOut[2]<0.0) {rgbOut[2]=0.0;} /* clip values up to zero */ \
\
          rgbOut[0] = gain * rgbOut[0]; \
          rgbOut[1] = gain * rgbOut[1]; \
          rgbOut[2] = gain * rgbOut[2]; \
\
          float rgb_value[3]; \
          rgb_value[0] = rgbOut[0]; \
          rgb_value[1] = rgbOut[1]; \
          rgb_value[2] = rgbOut[2]; \
\
/* compute simple fifth-power over fourth-power norm (yellow weighted norm might be less noise-prone, but stay simple for now) */ \
          float odt_norm = (rgbOut[0] * rgbOut[0] * rgbOut[0] * rgbOut[0] * rgbOut[0] + \
                            rgbOut[1] * rgbOut[1] * rgbOut[1] * rgbOut[1] * rgbOut[1] + \
                            rgbOut[2] * rgbOut[2] * rgbOut[2] * rgbOut[2] * rgbOut[2]) / \
                 MAX(TINY, (rgbOut[0] * rgbOut[0] * rgbOut[0] * rgbOut[0] + \
                            rgbOut[1] * rgbOut[1] * rgbOut[1] * rgbOut[1] + \
                            rgbOut[2] * rgbOut[2] * rgbOut[2] * rgbOut[2])); \
\
          float norm_asymp; \
\
/* make shoulder using asymptotic function */ \
\
          norm_asymp = odt_norm / (half_way + odt_norm); \
\
          norm_asymp = (1.0 - norm_asymp) *  norm_asymp /* low value behavior is just norm_asymp */ \
                           +  norm_asymp  * (norm_asymp + (1.0 - norm_asymp) * room_brightening); /* lift asymp for low-mids and mids */ \
\
          float asymp_scale = MAX(0.0, norm_asymp / (MAX(TINY, odt_norm))); \
\
          rgbOut[0] = asymp_scale * rgbOut[0]; \
          rgbOut[1] = asymp_scale * rgbOut[1]; \
          rgbOut[2] = asymp_scale * rgbOut[2]; \
\
/* allow gain to go higher, will reduce below after matrix, uncomment following lines (not print statements if OpenCL) to re-instate 1.0 norm ceiling */ \
/* rgbOut[0] = NORM_CORRECTION * rgbOut[0]; */ /* reduce gain based upon maximum norm value */ \
/* rgbOut[1] = NORM_CORRECTION * rgbOut[1]; */ /* reduce gain based upon maximum norm value */ \
/* rgbOut[2] = NORM_CORRECTION * rgbOut[2]; */ /* reduce gain based upon maximum norm value */ \
/* if (rgbOut[0]>1.0) {*/ /*printf("rgbOut[0]=%f, decrease NORM_CORRECTION a little\n", rgbOut[0]);*/ /*rgbOut[0] = 1.0; }*/ \
/* if (rgbOut[1]>1.0) {*/ /*printf("rgbOut[1]=%f, decrease NORM_CORRECTION a little\n", rgbOut[1]);*/ /*rgbOut[1] = 1.0; }*/ \
/* if (rgbOut[2]>1.0) {*/ /*printf("rgbOut[2]=%f, decrease NORM_CORRECTION a little\n", rgbOut[2]);*/ /*rgbOut[2] = 1.0; }*/ \
\
           rgb_value[0] = powf(rgb_value[0] / (half_way + rgb_value[0]), non_radiometric_asymp_gamma_boost); \
           rgb_value[1] = powf(rgb_value[1] / (half_way + rgb_value[1]), non_radiometric_asymp_gamma_boost); \
           rgb_value[2] = powf(rgb_value[2] / (half_way + rgb_value[2]), non_radiometric_asymp_gamma_boost); \
\
\
          float RGB_VEC[3], RGB_VECB[3]; \
\
          MATMUL(RGB_VEC, aces_to_device_mat, rgbOut) \
\
          MATMUL(RGB_VECB, aces_to_device_mat, rgb_value) \
\
/* blend radiometric chromaticity-preserving RGB_VEC and highlight-desaturationg RGB_VECB */ \
          RGB_VEC[0] = radiometric_proportion * RGB_VEC[0] + (1.0 - radiometric_proportion) * RGB_VECB[0]; \
          RGB_VEC[1] = radiometric_proportion * RGB_VEC[1] + (1.0 - radiometric_proportion) * RGB_VECB[1]; \
          RGB_VEC[2] = radiometric_proportion * RGB_VEC[2] + (1.0 - radiometric_proportion) * RGB_VECB[2]; \
\
/* clip values above 1.0, reducing the other primaries in relative proportion to the clip reduction to preserve xy chromaticity */ \
/* this > 1.0 clip will result in slope discontinuities.  It would be better to use a smooth function */ \
	  if(RGB_VEC[0] > 1.0) { RGB_VEC[1]=RGB_VEC[1]/RGB_VEC[0]; RGB_VEC[2]=RGB_VEC[2]/RGB_VEC[0]; RGB_VEC[0] = 1.0; } \
	  if(RGB_VEC[1] > 1.0) { RGB_VEC[0]=RGB_VEC[0]/RGB_VEC[1]; RGB_VEC[2]=RGB_VEC[2]/RGB_VEC[1]; RGB_VEC[1] = 1.0;} \
	  if(RGB_VEC[2] > 1.0) { RGB_VEC[0]=RGB_VEC[0]/RGB_VEC[2]; RGB_VEC[1]=RGB_VEC[1]/RGB_VEC[2]; RGB_VEC[2] = 1.0;} \
\
/* Clip negative RGB P3 Values, note that this is an out-of-gamut clip, meaning that xy chromaticity is not preserved */ \
/* this < 0 clip will also result in slope discontinuities.  It would be better to use a smooth function */ \
	  if(RGB_VEC[0] < 0) RGB_VEC[0] = 0; \
	  if(RGB_VEC[1] < 0) RGB_VEC[1] = 0; \
	  if(RGB_VEC[2] < 0) RGB_VEC[2] = 0; \
\
          process_gamma(display_gamma) \
}

/***********************************************************************************************************/

/* P3_D60_RGB_FROM_ACES */
#define LOAD_HDR_MATRIX \
          float HDR_MAT[9] = \
            { 1.980029295, -0.65352179,  -0.326511288, \
             -0.184520312,  1.288293303, -0.103777815, \
              0.008600852, -0.060073703,  1.05146064 };

/* P3_D65_RGB_FROM_ACES */
/* ACES_RGB_D60 to XYZ, then XYZ to P3_D65 */
/*#define LOAD_HDR_MATRIX \
          float HDR_MAT[9] = \
           { 2.054824895,   -.678201524,  -.338848384, \
             -.183835832,  1.283509252,  -.103389667, \
              .007944253,  -.055466159,   .970827554 }; */

/* ACES_RGB_D60_TO_BT2020_D60 */
/*#define LOAD_HDR_MATRIX \
          float HDR_MAT[9] = \
           { 1.47728,       -.252896,	   -.224381, \
             -.0792524,    1.18028,      -.101031, \
              .00226209,   -.0336893,     1.03143 }; */

/* ACES_RGB_D60_TO_BT2020_D65 */
/*#define LOAD_HDR_MATRIX \
          float HDR_MAT[9] = \
            { 1.512860179,   -.258987564,  -.229785733, \
              -.079036555,  1.177065932,  -.100755613, \
                .002091278,  -.031144101,   .953504037 };*/


#define process_odt_gd9_p3_d60_g2pt4_HDR \
{/* odt_type P3D60_HDR */ \
\
/* input is rgbOut[3], output is rOut, gOut, bOut */ \
\
   LOAD_HDR_MATRIX \
\
   xy_preserving_odt(HDR_GAIN, HALF_WAY_HDR, GAMMA_BOOST_HDR, ROOM_BRIGHTENING, HDR_MAT, HDR_DISPLAY_GAMMA, \
    HDR_RADIOMETRIC_PROPORTION) \
\
}

/***********************************************************************************************************/
/* REC709_RGB_D60_AS_D65_FROM_ACES */
#define LOAD_MDR_MATRIX \
          float MDR_MAT[9] = \
            { 2.5216353,  -1.136895218, -0.384904042, \
             -0.275205678, 1.369705042, -0.094399363, \
             -0.01593043, -0.147809266,  1.163803579 };


/* ACES_RGB_D60_TO_REC709_D65 (less magnitude in blue row) */
/*#define LOAD_MDR_MATRIX \
          float MDR_MAT[9] = \
            { 2.55839,      -1.11947,      -.391812, \
              -.277986,     1.36602,      -.0934873, \
              -.0171707,    -.148529,     1.08102 };*/


#define process_odt_gd9_Rec709_g2pt4_MDR \
{/* odt_type Rec709_MDR */ \
\
/* input is rgbOut[3] */ \
   LOAD_MDR_MATRIX \
\
   xy_preserving_odt(MDR_GAIN, HALF_WAY_HDR, GAMMA_BOOST_MDR, ROOM_BRIGHTENING, MDR_MAT, \
    MDR_DISPLAY_GAMMA, MDR_RADIOMETRIC_PROPORTION) \
\
}


#else /* not RADIOMETRIC_ODTS */
/***********************************************************************************************************/
#define process_odt_gd9_p3_d60_g2pt4_HDR \
{ /* odt_type P3D60_HDR */ \
\
/* boost the saturation a little */ \
\
          float SAT_MTX[9]; \
          float rgb_output[3]; \
          float off_diag_hdr_amount = HDR_GAIN *        OFF_DIAG_HDR; \
          float     diag_hdr_amount = HDR_GAIN * (1.0 - OFF_DIAG_HDR - OFF_DIAG_HDR); \
          SAT_MTX[0] = SAT_MTX[4] = SAT_MTX[8] = diag_hdr_amount; \
          SAT_MTX[1] = SAT_MTX[2] = SAT_MTX[3] = SAT_MTX[5] = SAT_MTX[6] = SAT_MTX[7] = off_diag_hdr_amount; \
\
          COLOR_MATRIX_TIMES_VEC(rgb_output, SAT_MTX, rgbOut) \
\
          rgbOut[0] = rgb_output[0]; /* can copy back once matrix transform is complete */ \
          rgbOut[1] = rgb_output[1]; \
          rgbOut[2] = rgb_output[2]; \
\
          float sign_of_rgbOut[3]; \
          for (i=0; i<3; i++) { \
            if (rgbOut[i] >= 0.0) { \
              sign_of_rgbOut[i] = 1.0; \
            } else { /* negative */ \
              sign_of_rgbOut[i] = -1.0; \
            } /* positive or not */ \
/* make shoulder using asymptotic function */ \
            rgbOut[i] = sign_of_rgbOut[i] * powf(sign_of_rgbOut[i] * rgbOut[i] / MAX(TINY, (HALF_WAY_HDR + (sign_of_rgbOut[i] * rgbOut[i]))), GAMMA_BOOST_HDR); \
          } /* i */ \
\
          float RGB_VEC[3]; \
/* ACES to XYZ */ \
/*        float XYZ_FROM_ACES_MAT[9], RGB_FROM_XYZ_MAT[9]; \
          XYZ_FROM_ACES_MAT[0] =  0.9525523959; \
          XYZ_FROM_ACES_MAT[1] =  0.0000000000; \
          XYZ_FROM_ACES_MAT[2] =  0.0000936786; \
          XYZ_FROM_ACES_MAT[3] =  0.3439664498; \
          XYZ_FROM_ACES_MAT[4] =  0.7281660966; \
          XYZ_FROM_ACES_MAT[5] = -0.0721325464; \
          XYZ_FROM_ACES_MAT[6] =  0.0000000000; \
          XYZ_FROM_ACES_MAT[7] =  0.0000000000; \
          XYZ_FROM_ACES_MAT[8] =  1.0088251844; \
\
          float XYZ_VEC[3]; \
          MATMUL(XYZ_VEC, XYZ_FROM_ACES_MAT, rgbOut) */\
\
/* matrix from XYZ to p3_D60 whitepoint */ \
/*          RGB_FROM_XYZ_MAT[0] = 2.40274; \
          RGB_FROM_XYZ_MAT[1] = -.89749; \
          RGB_FROM_XYZ_MAT[2] = -.38805; \
          RGB_FROM_XYZ_MAT[3] = -.83258; \
          RGB_FROM_XYZ_MAT[4] = 1.76923; \
          RGB_FROM_XYZ_MAT[5] =  .02371; \
          RGB_FROM_XYZ_MAT[6] =  .03882; \
          RGB_FROM_XYZ_MAT[7] = -.08250; \
          RGB_FROM_XYZ_MAT[8] = 1.03636; \
\
          MATMUL(RGB_VEC, RGB_FROM_XYZ_MAT, XYZ_VEC) */ \
\
/* net ACES to P3_D60 matrix by concatenating both matrices above: */ \
          float P3_D60_RGB_FROM_ACES_MAT[9]; \
          P3_D60_RGB_FROM_ACES_MAT[0] = 1.980029295; \
          P3_D60_RGB_FROM_ACES_MAT[1] = -0.65352179; \
          P3_D60_RGB_FROM_ACES_MAT[2] = -0.326511288; \
          P3_D60_RGB_FROM_ACES_MAT[3] = -0.184520312; \
          P3_D60_RGB_FROM_ACES_MAT[4] = 1.288293303; \
          P3_D60_RGB_FROM_ACES_MAT[5] = -0.103777815; \
          P3_D60_RGB_FROM_ACES_MAT[6] = 0.008600852; \
          P3_D60_RGB_FROM_ACES_MAT[7] = -0.060073703; \
          P3_D60_RGB_FROM_ACES_MAT[8] = 1.05146064; \
\
          MATMUL(RGB_VEC, P3_D60_RGB_FROM_ACES_MAT, rgbOut) \
\
/* Clip negative RGB P3 Values */ \
	  if(RGB_VEC[0] < 0) RGB_VEC[0] = 0; \
	  if(RGB_VEC[1] < 0) RGB_VEC[1] = 0; \
	  if(RGB_VEC[2] < 0) RGB_VEC[2] = 0; \
/* clip values above 1.0 */ \
	  if(RGB_VEC[0] > 1.0) RGB_VEC[0] = 1.0; \
	  if(RGB_VEC[1] > 1.0) RGB_VEC[1] = 1.0; \
	  if(RGB_VEC[2] > 1.0) RGB_VEC[2] = 1.0; \
\
          process_gamma(HDR_DISPLAY_GAMMA) \
}

/***********************************************************************************************************/

#define process_odt_gd9_Rec709_g2pt4_MDR \
{ /* odt_type Rec709_MDR */ \
\
/* boost the saturation a little */ \
\
          float SAT_MTX[9]; \
          float rgb_output[3]; \
          float off_diag_mdr_amount = MDR_GAIN *        OFF_DIAG_MDR; \
          float     diag_mdr_amount = MDR_GAIN * (1.0 - OFF_DIAG_MDR - OFF_DIAG_MDR); \
          SAT_MTX[0] = SAT_MTX[4] = SAT_MTX[8] = diag_mdr_amount; \
          SAT_MTX[1] = SAT_MTX[2] = SAT_MTX[3] = SAT_MTX[5] = SAT_MTX[6] = SAT_MTX[7] = off_diag_mdr_amount; \
\
          COLOR_MATRIX_TIMES_VEC(rgb_output, SAT_MTX, rgbOut) \
\
          rgbOut[0] = rgb_output[0]; /* can copy back once matrix transform is complete */ \
          rgbOut[1] = rgb_output[1]; \
          rgbOut[2] = rgb_output[2]; \
\
          float sign_of_rgbOut[3]; \
          for (i=0; i<3; i++) { \
            if (rgbOut[i] >= 0.0) { \
              sign_of_rgbOut[i] = 1.0; \
            } else { /* negative */ \
              sign_of_rgbOut[i] = -1.0; \
            } /* positive or not */ \
/* make shoulder using asymptotic function */ \
            rgbOut[i] = sign_of_rgbOut[i] * powf(sign_of_rgbOut[i] * rgbOut[i] / MAX(TINY, (HALF_WAY_MDR + (sign_of_rgbOut[i] * rgbOut[i]))), GAMMA_BOOST_MDR); \
          } /* i */ \
\
          float RGB_VEC[3]; \
\
/* ACES_RGB_D60 to Rec709_D60_as_D65: */ \
/* note that all of the rows of this matrix sum to 1.0, such that ACES=OCES=r=g=b=D60 is turned into Rec709 r=g=b=D65 */ \
\
          float RGB_FROM_ACES_MAT[9]; \
          RGB_FROM_ACES_MAT[0] =  2.5216353; \
          RGB_FROM_ACES_MAT[1] = -1.136895218; \
          RGB_FROM_ACES_MAT[2] = -0.384904042; \
          RGB_FROM_ACES_MAT[3] = -0.275205678; \
          RGB_FROM_ACES_MAT[4] =  1.369705042; \
          RGB_FROM_ACES_MAT[5] = -0.094399363; \
          RGB_FROM_ACES_MAT[6] = -0.01593043; \
          RGB_FROM_ACES_MAT[7] = -0.147809266; \
          RGB_FROM_ACES_MAT[8] =  1.163803579; \
          MATMUL(RGB_VEC, RGB_FROM_ACES_MAT, rgbOut) \
\
/* ACES_RGB_D60 to between P3_D60 and Rec709_D60_as_D65: */ \
/* note that all of the rows of this matrix sum to 1.0, such that ACES=OCES=r=g=b=D60 is turned into r=g=b=D65 */ \
/*          float RGB_FROM_ACES_MAT[9]; \
          RGB_FROM_ACES_MAT[0] =  2.25; \
          RGB_FROM_ACES_MAT[1] = -.95; \
          RGB_FROM_ACES_MAT[2] = -.3; \
          RGB_FROM_ACES_MAT[3] = -0.23; \
          RGB_FROM_ACES_MAT[4] =  1.33; \
          RGB_FROM_ACES_MAT[5] = -0.1; \
          RGB_FROM_ACES_MAT[6] = -0.01; \
          RGB_FROM_ACES_MAT[7] = -0.09; \
          RGB_FROM_ACES_MAT[8] =  1.10; \
          MATMUL(RGB_VEC, RGB_FROM_ACES_MAT, rgbOut) */\
\
/* Clip negative RGB Values */ \
	  if(RGB_VEC[0] < 0) RGB_VEC[0] = 0; \
	  if(RGB_VEC[1] < 0) RGB_VEC[1] = 0; \
	  if(RGB_VEC[2] < 0) RGB_VEC[2] = 0; \
/* clip values above 1.0 */ \
	  if(RGB_VEC[0] > 1.0) RGB_VEC[0] = 1.0; \
	  if(RGB_VEC[1] > 1.0) RGB_VEC[1] = 1.0; \
	  if(RGB_VEC[2] > 1.0) RGB_VEC[2] = 1.0; \
\
          process_gamma(MDR_DISPLAY_GAMMA) \
}

#endif /* RADIOMETRIC_ODTS or not */
/***********************************************************************************************************/

 #define ACES_ADJUST_LMT \
  { \
  float rr, gg, bb; \
  float sign_rt, sign_gt, sign_bt; \
  float red_pos, grn_pos, blu_pos; \
\
/* rotation matrix from ACES to P3 primary directions (ACES_RGB_D60 to P3_D60): */ \
      red_pos=aces[0]; grn_pos=aces[1]; blu_pos=aces[2]; \
      if (red_pos < 0.0) { red_pos =0.0; /* eliminate negative aces values being folded positive at negative matrix terms */ } \
      if (grn_pos < 0.0) { grn_pos =0.0; /* eliminate negative aces values being folded positive at negative matrix terms */ } \
      if (blu_pos < 0.0) { blu_pos =0.0; /* eliminate negative aces values being folded positive at negative matrix terms */ } \
\
      rr = 1.98003    * aces[0] -  .65352   * grn_pos -  .32651   * blu_pos; \
      gg = -.18452    * red_pos + 1.288298  * aces[1] -  .103778  * blu_pos; \
      bb =  .00860085 * aces[0] -  .0600737 * grn_pos + 1.0514732 * aces[2]; \
\
      rt = r_gain * rr; \
      gt = g_gain * gg; \
      bt = b_gain * bb; \
      rt = rt + r_offset; \
      gt = gt + g_offset; \
      bt = bt + b_offset; \
\
      sign_rt = 1.0; \
      if (rt < 0.0) { sign_rt = -1.0; } \
      sign_gt = 1.0; \
      if (gt < 0.0) { sign_gt = -1.0; } \
      sign_bt = 1.0; \
      if (bt < 0.0) { sign_bt = -1.0; } \
\
      rt = sign_rt * powf(sign_rt * rt, cdl_gamma_r); \
      gt = sign_gt * powf(sign_gt * gt, cdl_gamma_g); \
      bt = sign_bt * powf(sign_bt * bt, cdl_gamma_b); \
\
      if (saturation != 1.0) { /* process saturation using asc-cdl_1.2 formula */ \
/*	"sat" range: 0.0 <= sat < infinity \
	0.0 -> fully desaturated \
	1.0 -> no change \
	>1.0 -> more saturation */ \
\
/* compute the fully desaturated grey value */ \
        luma = rt * 0.2126 + gt * 0.7152 + bt * 0.0722; \
\
/*	output rgb is a weighted average between the input rgb and the fully desaturated grey value with "sat" as the weight... */ \
        rt = rt * saturation + luma * (1.0 - saturation); \
        gt = gt * saturation + luma * (1.0 - saturation); \
        bt = bt * saturation + luma * (1.0 - saturation); \
\
      } /* saturation != 1.0 */ \
\
/* rotation back from P3 direction primaries back to ACES primaries (P3_D60 to ACES_RGB_D60, inverse of above): */ \
      rr =     .53010 * rt +    .27786 * gt + .19204 * bt; \
      gg =     .07593 * rt +    .81961 * gt + .10446 * bt; \
      bb = /* 0.0     * rt + */ .04455 * gt + .95545 * bt; \
\
      aces[0] = rr; \
      aces[1] = gg; \
      aces[2] = bb; \
  } /* ACES_ADJUST_LMT */


/*************************************************************************************************/
#define compute_slow_in_out_360(colmin, colmax) \
        slow_in_out = 0.5 + 0.5 * sinf(2.0 * value_of_pi * (MAX(0.0, MIN(1.0, colmin/MAX(1e-9, colmax)))) - value_of_pi/2.0);

#define compute_slow_in_out_180(colmin, colmax) \
        slow_in_out = 0.5 + 0.5 * sinf(value_of_pi * (MAX(0.0, MIN(1.0, colmin/MAX(1e-9, colmax)))) + value_of_pi/2.0);


 #ifdef FACES_HIGHLIGHTS
  #define BRIGHT_MAX_FOR_SAT_BOOST .425//uuu.325
  #define BRIGHT_PWR_FOR_SAT_BOOST .425//uuu.5
  #define SAT_BOOST_AMOUNT .37//mmmm.34//mmm.45/*hhh.65*//*eee.5*//*eee0.415*/ /* note: this only applies to FULL_LMT */
  #define DESAT_BLEND_WEIGHT .75/*ooo.65*/ /* amount of desat sat boost vs amount of mid-sat sat boost, 1.0 is full desat (low and mid) sat boost, 0.0 is full mid_sat sat boost */
 #else /* not FACES_HIGHLIGHTS */
  #define BRIGHT_MAX_FOR_SAT_BOOST .975//uuu1.2
  #define BRIGHT_PWR_FOR_SAT_BOOST .36//uuu.33
  #define SAT_BOOST_AMOUNT .45//uuu0.35 /* note: this only applies to FULL_LMT */
  #define DESAT_BLEND_WEIGHT .75/*ooo.65*/ /* amount of desat sat boost vs amount of mid-sat sat boost, 1.0 is full desat (low and mid) sat boost, 0.0 is full mid_sat sat boost */
 #endif /* FACES_HIGHLIGHTS or not */

/***********************************************************************************************************/
// Algorithm for applying tone curve in forward direction.

/* modified rrt_shaper_fwd_gd10: */
/* this version by Gary Demos, Scott Dyer, Doug Walker */
/* this version slope-matched extended straight line in log-log space (i.e. a gamma) starting at BRIGHT_END, with BRIGHT_END being adjustable */
/* designed so that slope is linear .02 at low values (going down through zero to negative */
/* designed so that slope above BRIGHT_END tends toward unity slope (out = in) */

#define ACESMAX 65504.0
#define ACESMID 0.18
#define ACESMIN ACESMID * ACESMID / ACESMAX /* this ends up being .0000004946263 */
/* This is the same as: pow10(log10(ACESMID) - (log10(ACESMAX)-log10(ACESMID))); */
#define OCESMAX 10000.0
#define OCESMID 5.0
#define OCESMIN 0.0001

#define DARK_START -3.5 /* This is not easily adjustable.  COEFS 5 and 6 were matched ad-hoc to unity gamma for dark slope starting at DARK_START.  The match is slightly imperfect and discontinuous at present */
#define TINY_FOR_NORMS .000316228 /* matches DARK_START of -3.5 */
#define BRIGHT_END 1.5//2.5 /* adjust this to change the extended straight-line slope (i.e. a gamma) for the bright end extension, should end in .5 (half way between knots) */

#define RRT_KNOT_LEN 21

#define RRT_COEFS {-4.00000, -4.00000, -4.00000, -3.95000, -3.82000 /*<-ignored with dark_start at -3.5*/, -3.8, -3.245, -2.56, -1.8, -.825, .235,  1.05/*ppp2 1.075*//*ppp 1.125*/, 1.8 /*ppp2 1.925*/ /*ppp 1.94320*/, 2.55 /*ppp2 2.45*/ /*ppp 2.61950*/, 3.05 /*ppp4 2.87*//*ppp3 3.0*//*ppp2 3.195*/ /*ppp 3.13250*/, 3.6738 /*ppp 3.49050*/,  3.71550,  3.85130,  3.92710,  3.96980, 4.00000,  4.00000, 4.00000}

#define interp(value, out) \
    knot_coord = ( value - RRT_KNOT_START) / RRT_KNOT_SPAN * (RRT_KNOT_LEN-1); \
    j = knot_coord; \
    t = knot_coord - j; \
\
    wgt2 = 0.5 * coef[j] - coef[j+1] + 0.5 * coef[j+2]; \
    wgt1 = -coef[j] + coef[j+1]; \
    wgt0 = 0.5 * coef[j] + 0.5 * coef[j+1]; \
\
    out = t * t * wgt2 + t * wgt1 + wgt0;


/* Input and output are in linear (not log) units. */
#define rrt_shaper_fwd_gd10(aces_x, oces_y) \
{ \
\
  float knot_coord, t, wgt0, wgt1, wgt2; \
  short j; \
\
  float coef[23] = RRT_COEFS; \
\
  const float RRT_KNOT_START = log10f(ACESMIN); \
  const float RRT_KNOT_END = log10f(ACESMAX); \
  const float RRT_KNOT_SPAN = RRT_KNOT_END - RRT_KNOT_START; \
  const float RRT_KNOT_INC = RRT_KNOT_SPAN / (RRT_KNOT_LEN - 1.); \
\
  /* Check for negatives or zero before taking the log. If negative or zero, set to ACESMIN */ \
  float acesCheck = aces_x; \
  if (acesCheck <= ACESMIN) acesCheck = ACESMIN; \
\
  float logAces = log10f( acesCheck); \
\
  float logOces; \
\
  /* For logOces values in the knot range, apply the B-spline shaper, b(x) */ \
  if (( logAces >= DARK_START) && ( logAces < BRIGHT_END)) { \
    interp(logAces, logOces) \
    oces_y = powf(10.0, logOces); \
  } else if ( logAces < DARK_START) { \
    oces_y = aces_x; \
  } else if ( logAces >= BRIGHT_END) { \
\
\
    float log_end, log_end_fwd, log_end_back, log_slope; \
    float brt_end = BRIGHT_END; /* locations are half way in-between knots */ \
    float brt_end_fwd  = BRIGHT_END + 0.5 * RRT_KNOT_INC; \
    float brt_end_back = BRIGHT_END - 0.5 * RRT_KNOT_INC; \
\
/* note: the basic parameters here should be done in a setup call, such that they need not be recomputed every pixel */ \
    interp(brt_end_fwd, log_end_fwd) /* go forward half a knot */ \
    interp(brt_end,     log_end) /* on center */ \
    interp(brt_end_back, log_end_back) /* go back half a knot */ \
    log_slope = (log_end_fwd - log_end_back) / RRT_KNOT_INC; \
    logOces = log_end + (logAces - brt_end) * log_slope; /* project forward */ \
/*printf(" log_end_fwd = %f log_end = %f log_end_back = %f log_slope = %f RRT_KNOT_INC = %f logAces = %f logOces = %f\n", log_end_fwd, log_end, log_end_back, log_slope, RRT_KNOT_INC, logAces, logOces);*/ \
    oces_y = powf(10.0, logOces); \
  } \
\
    oces_y = oces_y * .02; /* scale down so that 5nits = 0.10 for LAD mid-grey */ \
\
} /* rrt_shaper_fwd_gd10 */


/************************************************************************************************************************************/
/* tilt the brightness toward yellow */
 #define RED_NORM_GAIN 1.22
 #define GRN_NORM_GAIN 1.2
 #define RED_NORM_MULTIPLE 1.75
 #define GRN_NORM_MULTIPLE 1.5

#define weighted_yellow_third_norm(aces_val, third_norm) \
/* (r^3+r^3+g^3+g^3+g^3+b^3)/(r^2+r^2+g^2+g^2+g^2+b^2) third order norm. */ \
/* the weighting presumes that green noise is lowest, red noise is next lowest, and blue noise is highest */ \
/*  This should be reasonably robust to negative numbers, although there might be small slope discontinuities when r, g, or b terms cross zero */ \
{ \
float rr, gg, bb; \
rr=aces_val[0] * RED_NORM_GAIN; \
gg=aces_val[1] * GRN_NORM_GAIN; \
bb=aces_val[2] * (3.0 - RED_NORM_GAIN - GRN_NORM_GAIN); /* r, g, and b gains should sum to 3.0, averaging unity for each so that neutral will map directly on the tone-curve lookup */ \
      third_norm = (RED_NORM_MULTIPLE * rr * rr * rr + \
                    GRN_NORM_MULTIPLE * gg * gg * gg + \
                                        bb * bb * bb)/ \
         (MAX(TINY, RED_NORM_MULTIPLE * rr * rr + \
                    GRN_NORM_MULTIPLE * gg * gg + \
                                        bb * bb )); \
}

/************************************************************************************************************************************/
#define weighted_third_norm(aces_val, third_norm) \
/* (r^3+r^3+g^3+g^3+g^3+b^3)/(r^2+r^2+g^2+g^2+g^2+b^2) third order norm. */ \
/* the weighting presumes that green noise is lowest, red noise is next lowest, and blue noise is highest */ \
/*  This should be reasonably robust to negative numbers, although there might be small slope discontinuities when r, g, or b terms cross zero */ \
      third_norm = (RED_NORM_MULTIPLE * aces_val[0] * aces_val[0] * aces_val[0] + \
                    GRN_NORM_MULTIPLE * aces_val[1] * aces_val[1] * aces_val[1] + \
                                        aces_val[2] * aces_val[2] * aces_val[2])/ \
         (MAX(TINY, RED_NORM_MULTIPLE * aces_val[0] * aces_val[0] + \
                    GRN_NORM_MULTIPLE * aces_val[1] * aces_val[1] + \
                                        aces_val[2] * aces_val[2] ));


#define BLUE_BALANCE_SAT 1.93//1.87//1.8 /* rebalance for different distances from neutral for min/max */
#define GREEN_BALANCE_SAT 1.55//1.53//1.575 /* rebalance for different distances from neutral for min/max */

/************************************************************************************************************************************/
/* boost_sat simulates the increased color saturation caused by the middle steep slope part of the tone S-Curve curve */
/* with slow_in_out from low to bright_max, and slow_out from desat to full sat, sat is boosted at low to mid saturation (no need to boost sat if already saturated) */
#define boost_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) \
{ \
\
  float slow_in_out, slow_in_out_b, slow_in_out_c, norm_aces; \
\
  weighted_third_norm(aces, norm_aces) \
\
  if (norm_aces > TINY_FOR_NORMS) { \
    float sat, max_val, min_val; \
    float grn_weighted = GREEN_BALANCE_SAT*aces[1]; \
    float blu_weighted = BLUE_BALANCE_SAT *aces[2]; \
    max_val = MAX(TINY_FOR_NORMS, MAX(aces[0], MAX(grn_weighted, blu_weighted))); \
    min_val = MAX(TINY_FOR_NORMS, MIN(aces[0], MIN(grn_weighted, blu_weighted))); \
    if ((max_val == aces[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted))) { /* use two-term norm to mitigate blue noise, also check that min is not TINY */ \
       float inv_grn = 1.0/grn_weighted; \
       float inv_blu = 1.0/blu_weighted; \
       float inv_weighted_two_term_norm = (inv_grn * inv_grn + inv_grn * inv_grn + inv_blu * inv_blu) / \
                                          (inv_grn           + inv_grn           + inv_blu); /* (g^2 + g^2 + b^2) / (g + g + b) */ \
       min_val = 1.0 / inv_weighted_two_term_norm; \
    } /* (max_val == aces[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted)) */ \
    sat = ( max_val - min_val) / max_val; /* note: this has slope discontinuities */ \
    compute_slow_in_out_180(sat, sat_range ) /* for low-and-mid-sat */ \
    slow_in_out_b = powf(slow_in_out, sat_pwr); \
    compute_slow_in_out_360(sat, sat_range ) /* boost mid sat */ \
    slow_in_out_c = slow_in_out; \
    slow_in_out_b = DESAT_BLEND_WEIGHT * slow_in_out_b + (1.0 - DESAT_BLEND_WEIGHT) * slow_in_out_c; \
    compute_slow_in_out_360(norm_aces, bright_max) /* for mid-bright */ \
    slow_in_out = powf(slow_in_out, bright_pwr); \
  } else { /* magnitude_aces <= TINY */ \
    slow_in_out_b = 0.0; \
    slow_in_out = 0.0; \
  } /* norm_aces > TINY or not */ \
  float adjust_sat; \
  adjust_sat = slow_in_out_b /* for not saturated */ * slow_in_out /* for mid_bright */; /* only apply to mid-bright and desat-to-mid-sat *GD* */ \
\
  float sat_increase = sat_boost * adjust_sat; \
  if (sat_increase > 0.0) { \
    new_aces[0] = aces[0] + sat_increase * (aces[0] - norm_aces); \
    new_aces[1] = aces[1] + sat_increase * (aces[1] - norm_aces); \
    new_aces[2] = aces[2] + sat_increase * (aces[2] - norm_aces); \
  } else { /* no boost */ \
    new_aces[0] = aces[0]; \
    new_aces[1] = aces[1]; \
    new_aces[2] = aces[2]; \
  } /* sat_increase > 1.0 */ \
\
} /* boost_sat */

/************************************************************************************************************************************/
#define boost_blue_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) \
{ \
\
  float slow_in_out, slow_in_out_b, norm_aces; \
\
  weighted_third_norm(aces, norm_aces) \
  if ((norm_aces > TINY_FOR_NORMS) && (aces[2] > aces[1]) && (aces[2] > aces[0])) { \
    float sat, max_val, min_val; \
    max_val = aces[2]; \
    min_val = MAX(TINY_FOR_NORMS, (aces[0]*aces[0]*aces[0] + aces[1]*aces[1]*aces[1])/(aces[0]*aces[0]+aces[1]*aces[1])); /* red and green norm */ \
    sat = ( max_val - min_val) / max_val;  \
    compute_slow_in_out_360(sat, sat_range ) /* for mid-sat */ \
    slow_in_out_b = powf(slow_in_out, sat_pwr); \
    compute_slow_in_out_180(norm_aces, bright_max) /* for mid-bright */ \
    slow_in_out = powf((1.0 - slow_in_out), bright_pwr); \
  } else { /* magnitude_aces <= TINY */ \
    slow_in_out_b = 0.0; \
    slow_in_out = 0.0; \
  } /* norm_aces > TINY or not */ \
  float adjust_sat; \
  adjust_sat = slow_in_out_b /* for not saturated */ * slow_in_out /* for mid_bright */; /* only apply to mid-bright and desat-to-mid-sat *GD* */ \
\
  float sat_increase = sat_boost * adjust_sat; \
  if (sat_increase > 0.0) { \
    new_aces[0] = aces[0] + sat_increase * (aces[0] - norm_aces); \
    new_aces[1] = aces[1] + sat_increase * (aces[1] - norm_aces); \
    new_aces[2] = aces[2] + sat_increase * (aces[2] - norm_aces); \
  } else { /* no boost */ \
    new_aces[0] = aces[0]; \
    new_aces[1] = aces[1]; \
    new_aces[2] = aces[2]; \
  } /* sat_increase > 1.0 */ \
\
} /* boost_blue_sat */

/************************************************************************************************************************************/
#define boost_blue_mid_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) \
{ \
\
  float slow_in_out, slow_in_out_b, norm_aces; \
\
  weighted_third_norm(aces, norm_aces) \
  if ((norm_aces > TINY_FOR_NORMS) && (aces[2] > aces[1]) && (aces[2] > aces[0])) { \
    float sat, max_val, min_val; \
    max_val = aces[2]; \
    min_val = MAX(TINY_FOR_NORMS, (aces[0]*aces[0]*aces[0] + aces[1]*aces[1]*aces[1])/(aces[0]*aces[0]+aces[1]*aces[1])); /* red and green norm */ \
    sat = ( max_val - min_val) / max_val;  \
    compute_slow_in_out_360(sat, sat_range ) /* for mid-sat */ \
    slow_in_out_b = powf(slow_in_out, sat_pwr); \
    compute_slow_in_out_360(norm_aces, bright_max) /* for mid-bright */ \
    slow_in_out = powf(slow_in_out, bright_pwr); \
  } else { /* magnitude_aces <= TINY */ \
    slow_in_out_b = 0.0; \
    slow_in_out = 0.0; \
  } /* norm_aces > TINY or not */ \
  float adjust_sat; \
  adjust_sat = slow_in_out_b /* for not saturated */ * slow_in_out /* for mid_bright */; /* only apply to mid-bright and desat-to-mid-sat *GD* */ \
\
  float sat_increase = sat_boost * adjust_sat; \
  if (sat_increase > 0.0) { \
    new_aces[0] = aces[0] + sat_increase * (aces[0] - norm_aces); \
    new_aces[1] = aces[1] + sat_increase * (aces[1] - norm_aces); \
    new_aces[2] = aces[2] + sat_increase * (aces[2] - norm_aces); \
  } else { /* no boost */ \
    new_aces[0] = aces[0]; \
    new_aces[1] = aces[1]; \
    new_aces[2] = aces[2]; \
  } /* sat_increase > 1.0 */ \
\
} /* boost_blue_mid_sat */

#define BLUE_BALANCE_FOR_DARKEN_SAT .9//.825 /* rebalance for different distances from neutral for min/max, larger targets blue and cyan */
#define GREEN_BALANCE_FOR_DARKEN_SAT .6//.9 /* rebalance for different distances from neutral for min/max */
#define YEL_THRESH 0.35 /* amount of yellow to fully activate yellow (yellow not darkened above this amount */
#define YELLOW_BOOST 1.1 /* boost yellow brightness */

/************************************************************************************************************************************/
#define darken_sat(aces_start, bright_amt, sat_thresh_low, sat_thresh_high, darken_sat_pwr, darken_amt, new_aces) \
{ \
    float sat, max_val, min_val, slow_in_out, slow_in_out_b, amt; \
    float grn_weighted = GREEN_BALANCE_FOR_DARKEN_SAT*aces_start[1]; \
    float blu_weighted = BLUE_BALANCE_FOR_DARKEN_SAT*aces_start[2]; \
    max_val = MAX(LARGER_TINY, MAX(aces_start[0], MAX(grn_weighted, blu_weighted))); /* skew to blue, since blue is far from white */ \
    min_val = MAX(TINY, MIN(aces_start[0], MIN(grn_weighted, blu_weighted))); /* skew to blue, since blue is far from white */  \
    if ((max_val == aces_start[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted))) { /* use two-term norm to mitigate blue noise, also check that min is not TINY */ \
       float inv_grn = 1.0/grn_weighted; \
       float inv_blu = 1.0/blu_weighted; \
       float inv_weighted_two_term_norm = (inv_grn * inv_grn + inv_grn * inv_grn + inv_blu * inv_blu) / \
                                          (inv_grn           + inv_grn           + inv_blu); /* (g^2 + g^2 + b^2) / (g + g + b) */ \
       min_val = 1.0 / inv_weighted_two_term_norm; \
    } /* (max_val == aces_start[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted)) */ \
    sat = (max_val - min_val) / max_val; /* note: this has slope discontinuities */ \
    compute_slow_in_out_180((sat - sat_thresh_low), (sat_thresh_high - sat_thresh_low) ) /* for high sat */ \
    slow_in_out_b = powf((1.0 - slow_in_out/* invert to high sat */), darken_sat_pwr /* values above 1.0 focus more on fully sat values, below 1.0 affect somewhat desat as well*/); \
\
    compute_slow_in_out_180(max_val, bright_amt) /* dark values are noisy in sat, ramp in with bright_amt using max_val */ \
    slow_in_out = 1.0 - slow_in_out; /* for brightness to ramp in */ \
    amt = 1.0 - darken_amt *  slow_in_out_b * slow_in_out; \
/* select yellow so that saturated yellow is not darkened */ \
    float max_unbiased = MAX(LARGER_TINY, MAX(aces_start[0], MAX(aces_start[1], aces_start[2]))); \
    float blue_pos = MAX(TINY, aces_start[2]); \
    float yellow_amt = (MAX(0.0, aces_start[1] - blue_pos) * MAX(0.0, aces_start[0] - blue_pos)) / max_unbiased; \
    compute_slow_in_out_180(yellow_amt, YEL_THRESH) \
    amt = amt * slow_in_out + YELLOW_BOOST * (1.0 - slow_in_out); \
    new_aces[0] = amt * aces_start[0]; \
    new_aces[1] = amt * aces_start[1]; \
    new_aces[2] = amt * aces_start[2]; \
} /* darken_sat */

#define BLUE_BALANCE_BRIGHTEN_DARK_SAT 1.0
#define GREEN_BALANCE_BRIGHTEN_DARK_SAT 1.0

/************************************************************************************************************************************/
#define brighten_dark_sat(aces, bright_amt, sat_thresh_low, sat_thresh_high, brighten_sat_pwr, brighten_magnitude, new_aces) \
{ \
    float sat, max_val, min_val, slow_in_out, slow_in_out_b, amt; \
    float grn_weighted = GREEN_BALANCE_BRIGHTEN_DARK_SAT*aces[1]; \
    float blu_weighted = BLUE_BALANCE_BRIGHTEN_DARK_SAT*aces[2]; \
    max_val = MAX(LARGER_TINY, MAX(aces[0], MAX(grn_weighted, blu_weighted))); /* skew to blue, since blue is far from white */ \
    min_val = MAX(TINY, MIN(aces[0], MIN(grn_weighted, blu_weighted))); /* skew to blue, since blue is far from white */  \
    if ((max_val == aces[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted))) { /* use two-term norm to mitigate blue noise, also check that min is not TINY */ \
       float inv_grn = 1.0/grn_weighted; \
       float inv_blu = 1.0/blu_weighted; \
       float inv_weighted_two_term_norm = (inv_grn * inv_grn + inv_grn * inv_grn + inv_blu * inv_blu) / \
                                          (inv_grn           + inv_grn           + inv_blu); /* (g^2 + g^2 + b^2) / (g + g + b) */ \
       min_val = 1.0 / inv_weighted_two_term_norm; \
    } /* (max_val == aces[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted)) */ \
    sat = (max_val - min_val) / max_val; /* note: this has slope discontinuities */ \
    compute_slow_in_out_180((sat - sat_thresh_low), (sat_thresh_high - sat_thresh_low) ) /* for high sat */ \
    slow_in_out_b = powf((1.0 - slow_in_out/* invert to high sat */), brighten_sat_pwr /* values above 1.0 focus more on fully sat values, below 1.0 affect somewhat desat as well*/); \
\
    compute_slow_in_out_360(max_val, bright_amt) /* dark values are noisy in sat, ramp in and out from 0 to bright_amt using max_val */ \
    amt = 1.0 + brighten_magnitude *  slow_in_out_b * slow_in_out; \
    new_aces[0] = amt * aces[0]; \
    new_aces[1] = amt * aces[1]; \
    new_aces[2] = amt * aces[2]; \
}


/************************************************************************************************************************************/
/* toward sat red */
#define R_RG  .07/*hhh.05*//*eee0.0*/
#define R_RB  .075/*hhh.055*//*eee.01*/
#define R_GR  0.0//jjj-.04/*-.04 makes reds go toward magenta*//*rrr.04*//*ooo-.04*///-.02
#define R_GB -.04/*rrr.04*//*ooo-.04*/
#define R_BR -.1425//jjj-.13875/*ooo-.135*///-.17
#define R_BG -.035/*ooo-.03*/
/* toward sat yellow */
#define Y_RG -.05/*eee-.1*/
#define Y_RB -.225/*eee-.3*/
#define Y_GR -0.1/*rrr.15*//*ooo-.1*///-.05
#define Y_GB -0.1/*rrr.1*//*ooo-.06*/
#define Y_BR .033125/*ooo.05*///0.0//.1
#define Y_BG -.01/*ooo.025*///0.0//.05
/* toward sat green */
#define G_RG  -.15//hhh-.1
#define G_RB  -.2//hhh-.3
#define G_GR -.05
#define G_GB -.06
#define G_BR  -.2
#define G_BG  -.1
/* toward sat cyan */
#define C_RG -.05
#define C_RB  .03
#define C_GR -.175/*aaa-.35*//*ooo-.2*//*ooo-.075*//*ooo-.03*/
#define C_GB .05/*ooo-.04*/
#define C_BR -.175/*aaa-.35*//*ooo-.2*//*-.07*/
#define C_BG .1/*ooo-.1*//*-.013*/
/* toward sat blue */
#define B_RG 0.0//-.06
#define B_RB 0.0/*aaa-.15*//*ooo-.06*///-.1//.05
#define B_GR 0.0/*aaa-.075*//*ooo0.0*/
#define B_GB 0.0/*aaa-.075*//*ooo-.03*/
#define B_BR -.3//hhh-.15/*aaa-.3*//*ooo-.225*///-.2//-.13
#define B_BG -.24//hhh-.12/*aaa-.24*//*ooo.1*//*ooo-.075*///-.05//-.03
/* toward sat magenta */
#define M_RG -.03//hhh-.06
#define M_RB  .05
#define M_GR -.025//hhh-.05
#define M_GB -.03//hhh-.06
#define M_BR -.065//hhh-.13
#define M_BG -.015//hhh-.03
/* neutral */
#define W_RG  0.07//kkk.15//hhh.1
#define W_RB  0.0375//kkk.075//jjj.1//hhh0.05
#define W_GR  0.04//kkk.08//jjj.05//hhh0.0
#define W_GB  0.0625//kkk.125//hhh.1//hhh0.0
#define W_BR -.04//jjj-.0625//jjj-.08//hhh-.15//hhh0.0
#define W_BG -.04//jjj-.0625//jjj-.08//hhh-.15//hhh0.0


/* caution, BRIGHTEN_Y and BRIGHTEN_R are noise prone by brightening yellow, leaving the appearance of darkened blue grains (because blue is not brightened) */ \
/* noise prone on blue grains when brightening */
#define BRIGHTEN_Y -.0875//mmmm-0.125//mmm.115//jjj.135//hhh.165//hhh.0975/*ooo.08*///0.03 /* brighten Yellow */

#define BRIGHTEN_R .0125//mmm-.025/*jjj.03*///kkk.0675/*ooo.055*///0.03 /* brighten Red */
#define DARKEN_G -.26//mmmm-.16//hhh-.225/*ooo-.15*/ /* darken Green */
#define DARKEN_C -.16//kkk-.2/*eee-.1*//*ooo-.05*//*ooo-.075*///-.15 /* darken Cyan */
#define DARKEN_B .065//kkk.05/*ooo-.1*//*ooo-.01*//*ooo-.0175*///-.035 /* darken or brighten Blue */

#define RED_CENTER .57//mmmm.45//mmm0.7//jjj0.25
#define RED_COMPACT_PWR .79//mmm.86/*jjj.835*///jjj.87 /* values below 1.0 compact red toward RED_CENTER */
#define RED_WEIGHT .45//mmm.175//jjj.3//jjj0.9
#define YELLOW_CENTER .55//mmmm.4//mmm.8//jjj0.35
#define YELLOW_COMPACT_PWR .8//mmmm.75//mmm.8 /* values below 1.0 compact yellow toward YELLOW_CENTER */
#define YELLOW_WEIGHT .35//mmm0.01//jjj.05/*jjj.1625*//*.15*//*jjj.1*/ /* dont need as much yellow compaction when brightening */
#define BLUE_CENTER 0.4
#define BLUE_COMPACT_PWR .7 /* values below 1.0 compact blue toward BLUE_CENTER */
#define BLUE_WEIGHT 1.0
#define DIALOUT_ASYMP_BRIGHTNESS_HALF_WAY 0.5 /* asymptotic brightness dialout halfway point */
#define DIALUP_ASYMP_DARK_HALF_WAY .1//jjj0.075 /* asymptotic dark-to-brighter activation halfway point */

#define H_BLUE_BALANCE_MATRIX_ROLLOFF 1.0 /* rebalance for different distances from neutral for min/max */
#define H_GREEN_BALANCE_MATRIX_ROLLOFF 1.0 /* rebalance for different distances from neutral for min/max */

/************************************************************************************************************************************/
#define hue_dynamic_matrix(aces_source, bright_amt, sat_thresh_low_hue_mats, sat_thresh_high_hue_mats, sat_pwr_hue_mats, sat_thresh_low_dialout, sat_thresh_high_dialout, sat_pwr_dialout, newer_aces) \
{ \
    float recip_sum, delta_hue, sat, max_val, mid_val, min_val, slow_in_out, slow_in_out_b, slow_in_out_c, amt, rg_processed, rb_processed; \
    float wgt_r, wgt_y, wgt_g, wgt_c, wgt_b, wgt_m; \
    float red_weighted = MAX(TINY,                                aces_source[0]); \
    float grn_weighted = MAX(TINY, H_GREEN_BALANCE_MATRIX_ROLLOFF*aces_source[1]); \
    float blu_weighted = MAX(TINY, H_BLUE_BALANCE_MATRIX_ROLLOFF *aces_source[2]); \
    max_val = MAX(red_weighted, MAX(grn_weighted, blu_weighted)); \
    min_val = MIN(red_weighted, MIN(grn_weighted, blu_weighted));  \
    if ((max_val == red_weighted) && ((min_val == blu_weighted) || (min_val == grn_weighted)) && ((max_val - min_val) > TINY)) { /* red max, use two-term norm to mitigate blue noise, also check that max is at least TINY bigger than min (for division denominator) */ \
      if (min_val == blu_weighted) { \
/* r > g > b */ \
        mid_val = grn_weighted; \
        delta_hue =(mid_val - min_val) / (max_val - min_val); \
        compute_slow_in_out_180(delta_hue, 1.0) \
        wgt_r = slow_in_out; \
        wgt_y = 1.0 - slow_in_out; \
        wgt_g = 0.0; \
        wgt_c = 0.0; \
        wgt_b = 0.0; \
        wgt_m = 0.0; \
      } else { /* min_val == grn_weighted */ \
/* r > b > g */ \
        mid_val = blu_weighted; \
        delta_hue = (mid_val - min_val) / (max_val - min_val); \
        compute_slow_in_out_180(delta_hue, 1.0) \
        wgt_r = slow_in_out; \
        wgt_y = 0.0; \
        wgt_g = 0.0; \
        wgt_c = 0.0; \
        wgt_b = 0.0; \
        wgt_m = 1.0 - slow_in_out; \
      } /* min_val == blu_weighted or not */ \
      float inv_grn = 1.0/grn_weighted; \
      float inv_blu = 1.0/blu_weighted; \
      float inv_weighted_two_term_norm = (inv_grn * inv_grn + inv_grn * inv_grn + inv_blu * inv_blu) / \
                                         (inv_grn           + inv_grn           + inv_blu); /* (g^2 + g^2 + b^2) / (g + g + b) */ \
      min_val = 1.0 / inv_weighted_two_term_norm; \
    } else { /* red not max */ \
      if ((max_val == aces[2]) && ((min_val == red_weighted) || (min_val == grn_weighted)) && ((max_val - min_val) > TINY)) { /* blue max */ \
        if (min_val == red_weighted) { \
/* b > g > r */ \
          mid_val = grn_weighted; \
          delta_hue = (mid_val - min_val) / (max_val - min_val); \
          compute_slow_in_out_180(delta_hue, 1.0) \
          wgt_r = 0.0; \
          wgt_y = 0.0; \
          wgt_g = 0.0; \
          wgt_c = 1.0 - slow_in_out; \
          wgt_b = slow_in_out; \
          wgt_m = 0.0; \
        } else { /* min_val == grn_weighted */ \
/* b > r > g */ \
          mid_val = red_weighted; \
          delta_hue = (mid_val - min_val) / (max_val - min_val); \
          compute_slow_in_out_180(delta_hue, 1.0) \
          wgt_r = 0.0; \
          wgt_y = 0.0; \
          wgt_g = 0.0; \
          wgt_c = 0.0; \
          wgt_b = slow_in_out; \
          wgt_m = 1.0 - slow_in_out; \
        } /* min_val == red_weighted or not */ \
        float inv_red = 1.0/red_weighted; \
        float inv_grn = 1.0/grn_weighted; \
        float inv_weighted_two_term_norm = (inv_red * inv_red + inv_grn * inv_grn + inv_grn * inv_grn) / \
                                           (inv_red           + inv_grn           + inv_grn); /* (r^2 + g^2 + g^2) / (r + g + g) */ \
        min_val = 1.0 / inv_weighted_two_term_norm; \
      } else { /* blue not max */ \
        if ((max_val == grn_weighted) && ((min_val == red_weighted) || (min_val == blu_weighted)) && ((max_val - min_val) > TINY)) { /* grn max */ \
          if (min_val == red_weighted) { \
/* g > b > r */ \
            mid_val = blu_weighted; \
            delta_hue = (mid_val - min_val) / (max_val - min_val); \
            compute_slow_in_out_180(delta_hue, 1.0) \
            wgt_r = 0.0; \
            wgt_y = 0.0; \
            wgt_g = slow_in_out; \
            wgt_c = 1.0 - slow_in_out; \
            wgt_b = 0.0; \
            wgt_m = 0.0; \
          } else { /* min_val == blu_weighted */ \
/* g > r > b */ \
            mid_val = red_weighted; \
            delta_hue = (mid_val - min_val) / (max_val - min_val); \
            compute_slow_in_out_180(delta_hue, 1.0) \
            wgt_r = 0.0; \
            wgt_y = 1.0 - slow_in_out; \
            wgt_g = slow_in_out; \
            wgt_c = 0.0; \
            wgt_b = 0.0; \
            wgt_m = 0.0; \
          } /* min_val == red_weighted or not */ \
          float inv_red = 1.0/red_weighted; \
          float inv_blu = 1.0/blu_weighted; \
          float inv_weighted_two_term_norm = (inv_red * inv_red + inv_red * inv_red + inv_blu * inv_blu) / \
                                             (inv_red           + inv_red           + inv_blu); /* (r^2 + r^2 + b^2) / (r + r + b) */ \
          min_val = 1.0 / inv_weighted_two_term_norm; \
        } else { /* grn not max */ \
          wgt_r = 0.0; \
          wgt_y = 0.0; \
          wgt_g = 0.0; \
          wgt_c = 0.0; \
          wgt_b = 0.0; \
          wgt_m = 0.0; \
        } /* grn max or not */ \
      } /* blue max or not */ \
    } /* red max or not */ \
\
    sat = (max_val - min_val) / max_val; /* note: this has slope discontinuities */ \
    compute_slow_in_out_180((sat - sat_thresh_low_dialout), (sat_thresh_high_dialout - sat_thresh_low_dialout) ) /* for low sat */ \
    slow_in_out_b = powf(slow_in_out, sat_pwr_dialout /* values below 1.0 increase the influence of off-diagonal terms over the sat slow_out span */); \
\
    compute_slow_in_out_180((sat - sat_thresh_low_hue_mats), (sat_thresh_high_hue_mats - sat_thresh_low_hue_mats) ) /* for low sat */ \
    slow_in_out_c = powf(slow_in_out, sat_pwr_hue_mats /* values below 1.0 increase the influence of off-diagonal terms over the sat slow_out span */); \
\
    wgt_r = (1.0 - slow_in_out_c) * wgt_r; /* ramp up with increasing sat */ \
    wgt_y = (1.0 - slow_in_out_c) * wgt_y; /* ramp up with increasing sat */ \
    wgt_g = (1.0 - slow_in_out_c) * wgt_g; /* ramp up with increasing sat */ \
    wgt_c = (1.0 - slow_in_out_c) * wgt_c; /* ramp up with increasing sat */ \
    wgt_b = (1.0 - slow_in_out_c) * wgt_b; /* ramp up with increasing sat */ \
    wgt_m = (1.0 - slow_in_out_c) * wgt_m; /* ramp up with increasing sat */ \
\
    compute_slow_in_out_180(max_val, bright_amt) /* dark values are noisy in sat, ramp in from 0 to bright_amt using max_val */ \
    slow_in_out = 1.0 - slow_in_out; /* make bright values fully active, make dark values ramp off */ \
    amt = slow_in_out_b * slow_in_out; \
\
    float rg, rb, gr, gb, br, bg, ww; \
    ww = (1.0 - wgt_r - wgt_g - wgt_b); \
    rg = amt * (ww * W_RG + wgt_r * R_RG + wgt_y * Y_RG + wgt_g * G_RG + wgt_c * C_RG + wgt_b * B_RG + wgt_m * M_RG); \
    rb = amt * (ww * W_RB + wgt_r * R_RB + wgt_y * Y_RB + wgt_g * G_RB + wgt_c * C_RB + wgt_b * B_RB + wgt_m * M_RB); \
    gr = amt * (ww * W_GR + wgt_r * R_GR + wgt_y * Y_GR + wgt_g * G_GR + wgt_c * C_GR + wgt_b * B_GR + wgt_m * M_GR); \
    gb = amt * (ww * W_GB + wgt_r * R_GB + wgt_y * Y_GB + wgt_g * G_GB + wgt_c * C_GB + wgt_b * B_GB + wgt_m * M_GB); \
    br = amt * (ww * W_BR + wgt_r * R_BR + wgt_y * Y_BR + wgt_g * G_BR + wgt_c * C_BR + wgt_b * B_BR + wgt_m * M_BR); \
    bg = amt * (ww * W_BG + wgt_r * R_BG + wgt_y * Y_BG + wgt_g * G_BG + wgt_c * C_BG + wgt_b * B_BG + wgt_m * M_BG); \
\
\
newer_aces[0]= (1.0 - rg - rb) * aces_source[0] + \
                   rg * MAX(0.0, aces_source[1]/*discontinuous slope*/) + \
                   rb * MAX(0.0, aces_source[2]/*discontinuous slope*/); \
newer_aces[1]=     gr * MAX(0.0, aces_source[0]/*discontinuous slope*/) + \
               (1.0 - gr - gb) * aces_source[1] + \
                   gb * MAX(0.0, aces_source[2]/*discontinuous slope*/); \
newer_aces[2]=     br * MAX(0.0, aces_source[0]/*discontinuous slope*/) + \
                   bg * MAX(0.0, aces_source[1]/*discontinuous slope*/) + \
               (1.0 - br - bg) * aces_source[2]; \
\
/*wgt_y=0; \
wgt_r=0; \
wgt_g=0; \
wgt_c=0; \
wgt_b=0;*/ \
/* caution, wgt_y and wgt_r are noise prone, for example brightening yellow (which is anti-blue) gives the appearance of darkening blue grains (since blue is not brightened) */ \
  float bright_scaling = (1.0 + BRIGHTEN_Y * wgt_y + BRIGHTEN_R * wgt_r + DARKEN_G * wgt_g + DARKEN_C * wgt_c + DARKEN_B * wgt_b); \
  newer_aces[0] = newer_aces[0] * bright_scaling; \
  newer_aces[1] = newer_aces[1] * bright_scaling; \
  newer_aces[2] = newer_aces[2] * bright_scaling; \
/* compact red and blue brightness dynamic ranges */ \
  float bright_dialout_asymp = 1.0 - max_val / (max_val + DIALOUT_ASYMP_BRIGHTNESS_HALF_WAY); \
  float dark_dialout_asymp = max_val / (max_val + DIALUP_ASYMP_DARK_HALF_WAY); \
  float wgt_rr = RED_WEIGHT    * wgt_r * bright_dialout_asymp * dark_dialout_asymp; \
  float wgt_yy = YELLOW_WEIGHT * wgt_y * bright_dialout_asymp * dark_dialout_asymp; \
  float wgt_bb = BLUE_WEIGHT   * wgt_b * bright_dialout_asymp * dark_dialout_asymp; \
  newer_aces[0] = (1.0 - wgt_rr) * newer_aces[0] + wgt_rr * RED_CENTER     * powf(MAX(0.0, (newer_aces[0] / RED_CENTER)),     RED_COMPACT_PWR); \
  newer_aces[0] = (1.0 - wgt_yy) * newer_aces[0] + wgt_yy * YELLOW_CENTER  * powf(MAX(0.0, (newer_aces[0] / YELLOW_CENTER )), YELLOW_COMPACT_PWR); \
  newer_aces[1] = (1.0 - wgt_yy) * newer_aces[1] + wgt_yy * YELLOW_CENTER  * powf(MAX(0.0, (newer_aces[1] / YELLOW_CENTER )), YELLOW_COMPACT_PWR); \
  newer_aces[2] = (1.0 - wgt_bb) * newer_aces[2] + wgt_bb * BLUE_CENTER    * powf(MAX(0.0, (newer_aces[2] / BLUE_CENTER)),    BLUE_COMPACT_PWR); \
}

/************************************************************************************************************************************/
#define RG 0.007//jjj0.0//jjj-.10/*ttt-.12*//*ttt-.08*//*ttt-.04*///-.05//-.0375
#define RB -.025/*ttt-.05*//*ttt.02*///-.03//-.024
#define GREEN_IN_RED .075//iii.05//iii0.0/*hhh0.05*//*ttt-.05*///-.04//-.01//-.02//-.04//.01//-.01//-.03//.03
#define GREEN_IN_BLUE 0.06/*ttt-.06*///-.08//-.03 // if negative, this functions as blue reduces green */
#define RED_REDUCES_BLUE -.2//jjj-.15//jjj-.135/*ttt-.13*///-.12// /* adjusts faces yellow (more negative values) vs faces pink (less negative values), this is sensitive */
#define GREEN_REDUCES_BLUE -.1//jjj-.05//jjj-.03//-.03 /* if positive, this functions as blue in green */

#define BLUE_BALANCE_MATRIX_ROLLOFF 1.4//iii1.5 /* rebalance for different distances from neutral for min/max */
#define GREEN_BALANCE_MATRIX_ROLLOFF 1.25//iii1.2 /* rebalance for different distances from neutral for min/max */

#define lmt_dynamic_matrix(aces_source, bright_amt_low, bright_amt_high, sat_thresh_low, sat_thresh_high, sat_pwr, newer_aces) \
{ \
    float sat, max_val, min_val, slow_in_out, slow_in_out_b,  slow_in_out_c, amt; \
    float red_weighted = MAX(TINY, aces_source[0]); \
    float grn_weighted = MAX(TINY, GREEN_BALANCE_MATRIX_ROLLOFF*aces_source[1]); \
    float blu_weighted = MAX(TINY, BLUE_BALANCE_MATRIX_ROLLOFF*aces_source[2]); \
    max_val = MAX(red_weighted, MAX(grn_weighted, blu_weighted)); /* skew to blue, since blue is far from white */ \
    min_val = MIN(red_weighted, MIN(grn_weighted, blu_weighted)); /* skew to blue, since blue is far from white */  \
    if ((max_val == red_weighted) && ((min_val == blu_weighted) || (min_val == grn_weighted)) && ((max_val - min_val) > TINY)) { /* use two-term norm to mitigate blue noise, also check that min is not TINY */ \
       float inv_grn = 1.0/grn_weighted; \
       float inv_blu = 1.0/blu_weighted; \
       float inv_weighted_two_term_norm = (inv_grn * inv_grn + inv_grn * inv_grn + inv_blu * inv_blu) / \
                                          (inv_grn           + inv_grn           + inv_blu); /* (g^2 + g^2 + b^2) / (g + g + b) */ \
       min_val = 1.0 / inv_weighted_two_term_norm; \
    } /* (max_val == aces[0]) && ((min_val == blu_weighted) || (min_val == grn_weighted)) */ \
    sat = (max_val - min_val) / max_val; /* note: this has slope discontinuities */ \
    compute_slow_in_out_180((sat - sat_thresh_low), (sat_thresh_high - sat_thresh_low) ) /* for low sat */ \
    slow_in_out_b = powf(slow_in_out, sat_pwr /* values below 1.0 increase the influence of off-diagonal terms over the sat slow_out span */); \
\
    compute_slow_in_out_180(max_val, bright_amt_low) /* dark values are noisy in sat, ramp in from 0 to bright_amt using max_val */ \
    slow_in_out_c = 1.0 - slow_in_out; /* make bright values fully active, make dark values ramp off */ \
    compute_slow_in_out_180(max_val, bright_amt_high) /* ramp out matrices at bright values */ \
    amt = slow_in_out * slow_in_out_b * slow_in_out_c; \
float grn_in_red, grn_in_blu, red_reduces_blu, grn_reduces_blu, rg, rb; \
    rg =              amt * RG; \
    rb =              amt * RB; \
    grn_in_red =      amt * GREEN_IN_RED; \
    grn_in_blu =      amt * GREEN_IN_BLUE; \
    red_reduces_blu = amt * RED_REDUCES_BLUE; \
    grn_reduces_blu = amt * GREEN_REDUCES_BLUE; \
\
newer_aces[0]= (1.0 - rg - rb) * aces_source[0] + rg * MAX(0.0, aces_source[1]/*discontinuous, use if RG is negative */) + \
                rb * MAX(0.0, aces_source[2]/*discontinuous, use if RB is negative */); \
newer_aces[1]= grn_in_red * MAX(0.0, aces_source[0]/*discontinuous, use if GREEN_IN_RED is negative */) + \
             (1.0 - grn_in_red - grn_in_blu) * aces_source[1] + grn_in_blu * MAX(0.0, aces_source[2]/*discontinuous, use if GREEN_IN_BLUE is negative */); \
newer_aces[2]= red_reduces_blu * MAX(0.0, aces_source[0]/*discontinuous*/) + \
               grn_reduces_blu * MAX(0.0, aces_source[1]/*discontinuous, use if GREEN_REDUCES_BLUE is negative */) + \
               (1.0 - red_reduces_blu - grn_reduces_blu) * aces_source[2]; \
}


/************************************************************************************************************/
#define simple_processing(aces_in, aces_simple) \
  { \
/* Adjust ACES values */ \
      float aces_f[3], aces_e[3], aces_d[3], aces_c[3], aces_b[3], aces_a[3]; \
\
/* #define lmt_dynamic_matrix_for(aces_in, bright_amt, sat_thresh_low, sat_thresh_high, sat_pwr, aces_a) */ \
      lmt_dynamic_matrix(aces_in, .1, 1.5 * BRIGHT_MAX_FOR_SAT_BOOST, .625, 1.0, .65, aces_a) \
\
/*#define boost_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) */ \
      boost_sat(aces_a, 1.7, .6, BRIGHT_PWR_FOR_SAT_BOOST, BRIGHT_MAX_FOR_SAT_BOOST, .26/*ttt.325*//*.275*//*.2*/, aces_f) /* boost saturation for middle-bright */ \
      boost_sat(aces_f, 1.4, .7, .85, .2, .45, aces_simple) /* boost saturation in shadows */ \
  }
/************************************************************************************************************/
  #define RED_IN_GRN_GAMMA 0.1
  #define RED_FROM_BLU_GAMMA -.035
  #define GRN_FROM_BLU_GAMMA -.025

  #ifdef BRIGHT_HIGHLIGHTS_IN_TONE_CURVE
   #define GAMMA_BOOST_MAX 3.0 /* maximum norm value for gamma boost, at which point the reduced gamma is used */
  #else /* not BRIGHT_HIGHLIGHTS_IN_TONE_CURVE */
   #define GAMMA_BOOST_MAX 1.0 /* maximum norm value for gamma boost, at which point the reduced gamma is used */
  #endif /* BRIGHT_HIGHLIGHTS_IN_TONE_CURVE or not */

  #define GAMMA_FOR_GAM_AND_MAT_R 1.4//kkk1.35//kkk1.25
  #define GAMMA_FOR_GAM_AND_MAT_G 1.4//kkk1.35//kkk1.25
  #define GAMMA_FOR_GAM_AND_MAT_B 1.5//kkk1.45//kkk1.35//kkk1.25

  #define GAMMA_REDUCED_R 0.85//kkk0.8 /* reduced gamma for bright regions (for face highlights) */
  #define GAMMA_REDUCED_G 0.85//kkk0.8 /* reduced gamma for bright regions (for face highlights) */
  #define GAMMA_REDUCED_B 0.85//kkk0.8 /* reduced gamma for bright regions (for face highlights) */

  #define NORM_PWR_FOR_GAM_AND_MAT 0.9375//kkk0.95 /* values below 1.0 raise brightness of dark values, decrease brightness of bright values, beyond the tone curve */ 
  #define SCALE_FOR_GAM_AND_MAT .875//kkk.9

#define gamma_plus_mat(aces_in, aces_gamma_mat) \
{ float rr, sign_rr, gg, sign_gg, bb, sign_bb; \
   if (aces_in[0] < 0) { \
     sign_rr = -1.0; \
     rr = - aces_in[0]; \
   } else { \
     sign_rr = 1.0; \
     rr = aces[0]; \
   } \
   gg = (1.0 - RED_IN_GRN_GAMMA) * aces[1] + RED_IN_GRN_GAMMA * aces[0]; \
   if (gg < 0) { \
     sign_gg = -1.0; \
     gg = - gg; \
   } else { \
     sign_gg = 1.0; \
     /*gg = gg;*/ \
   } \
   bb = (1.0 - (RED_FROM_BLU_GAMMA) - (GRN_FROM_BLU_GAMMA)) * aces[2] + (RED_FROM_BLU_GAMMA) * MAX(0.0, aces[0]) + (GRN_FROM_BLU_GAMMA) * MAX(0.0, aces[1]); \
   if (bb < 0) { \
     sign_bb = -1.0; \
     bb = - bb; \
   } else { \
     sign_bb = 1.0; \
     /*bb = bb;*/ \
   } \
\
   float rgb_norm, gamma_r, gamma_g, gamma_b, slow_in_out; \
\
   rgb_norm = (rr*rr*rr + gg*gg*gg + bb*bb*bb)/MAX(TINY, (rr*rr + gg*gg + bb*bb)); \
\
   rgb_norm = MAX(TINY_FOR_NORMS, rgb_norm); \
\
   float rgb_norm_compacted; \
   rgb_norm_compacted = SCALE_FOR_GAM_AND_MAT * powf(rgb_norm, NORM_PWR_FOR_GAM_AND_MAT); /* raise brightness of dark values, decrease brightness of bright values, beyond the tone curve */ \
\
   compute_slow_in_out_180(rgb_norm, GAMMA_BOOST_MAX)  \
   gamma_r = (1.0 - slow_in_out) * GAMMA_REDUCED_R + slow_in_out * GAMMA_FOR_GAM_AND_MAT_R; /* fade out gamma_boost (down to 1.0) for bright norm values */ \
   gamma_g = (1.0 - slow_in_out) * GAMMA_REDUCED_G + slow_in_out * GAMMA_FOR_GAM_AND_MAT_G; /* fade out gamma_boost (down to 1.0) for bright norm values */ \
   gamma_b = (1.0 - slow_in_out) * GAMMA_REDUCED_B + slow_in_out * GAMMA_FOR_GAM_AND_MAT_B; /* fade out gamma_boost (down to 1.0) for bright norm values */ \
\
   aces_gamma_mat[0] = sign_rr * rgb_norm_compacted * powf(rr / rgb_norm, gamma_r); \
   aces_gamma_mat[1] = sign_gg * rgb_norm_compacted * powf(gg / rgb_norm, gamma_g); \
   aces_gamma_mat[2] = sign_bb * rgb_norm_compacted * powf(bb / rgb_norm, gamma_b); \
  }


/************************************************************************************************************/
#define moderate_processing(aces_in, aces_moderate) \
   { \
      float aces_f[3], aces_e[3], aces_d[3], aces_c[3], aces_b[3], aces_a[3]; \
\
/* #define lmt_dynamic_matrix_for(aces_in, bright_amt, sat_thresh_low, sat_thresh_high, sat_pwr, aces_a) */ \
      lmt_dynamic_matrix(aces_in, .1, 1.5 * BRIGHT_MAX_FOR_SAT_BOOST, .45/*iii.625*/, .85/*iii1.0*/, .65, aces_a) \
\
/*#define boost_blue_mid_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) */ \
      boost_blue_mid_sat(aces_a,  .7, .9, .7,  .75, .26/*jjj.35*//*ppp.3*/, aces_b) /* boost saturation of mid-bright blue */ \
/*#define darken_sat(aces, bright_amt, sat_thresh_low, sat_thresh_high, darken_sat_pwr, darken_amt, new_aces)*/ \
      darken_sat(aces_b, .85, .425, .95, 1.0, .15/*ppp.08*/, aces_c) /* darken saturated color brightness for mid-bright */ \
      darken_sat(aces_c, 2.25, .55, 1.0, 1.0, .44/*iii.35*/, aces_d) /* darken saturated color brightness for bright */ \
/*#define brighten_dark_sat(aces, bright_amt, sat_thresh_low, sat_thresh_high, brighten_sat_pwr, brighten_magnitude, new_aces)*/ \
      brighten_dark_sat(aces_d, .2/*ooo.25*/, 0.0/*ooo.275*/, .85/*ooo1.0*/, .8/*ooo.9*/, .125/*ooo.3*//*.24*//*.225*/,  aces_e) /* brighten dark saturated colors */ \
/*#define boost_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) */ \
      boost_sat(aces_e, 1.7, .6, BRIGHT_PWR_FOR_SAT_BOOST, BRIGHT_MAX_FOR_SAT_BOOST, .3/*iii.05*//*hhh.285*//*hhh.325*//*.275*//*.2*/, aces_f) /* boost saturation for middle-bright */ \
      boost_sat(aces_f, 1.4, .7, .85, .2, .3/*hhh.35*//*hhh.55*//*.45*/, aces_moderate) /* boost saturation in shadows */ \
   }


/************************************************************************************************************/
#define full_processing(aces_in, aces_full) \
   { \
      float aces_f[3], aces_e[3], aces_d[3], aces_c[3], aces_b[3], aces_aa[3], aces_a[3]; \
\
/* #define hue_dynamic_matrix_for(aces_in, bright_amt, sat_thresh_low_hue_mats, sat_thresh_high_hue_mats, sat_pwr_hue_mats, sat_thresh_low_dialout, sat_thresh_high_dialout, sat_pwr_dialout, aces_a) */ \
      hue_dynamic_matrix(aces_in, .1, .025/*jjj.15*//*.1*/, .45/*.4*/, 0.5/*jjj1.5*/, .43/*jjj.275*//*.25*/, 1.0, 1.0, aces_a) \
\
/*#define boost_blue_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) */ \
      boost_blue_sat(aces_a,     .75, 1.0, 1.0,  .75, .05/*aaa.12*/, aces_aa) /* boost saturation of blue */ \
      boost_blue_mid_sat(aces_aa, .7, 1.0,  .7,  .75, .055/*aaa.13*/, aces_b) /* boost saturation of mid-bright blue */ \
/*#define darken_sat(aces, bright_amt, sat_thresh_low, sat_thresh_high, darken_sat_pwr, darken_amt, new_aces)*/ \
      darken_sat(aces_b, .85, .425, .95, 1.0, .35/*mmmm.2*//*mmm.04*//*jjj.115*//*ooo.15*//*ppp.08*/, aces_c) /* darken saturated color brightness for mid-bright */ \
      darken_sat(aces_c, 2.25, .55, 1.0, 1.0, .35/*mmmm.2*//*mmm.075*//*jjj.263*//*ooo.35*/, aces_d) /* darken saturated color brightness for bright */ \
/*#define brighten_dark_sat(aces, bright_amt, sat_thresh_low, sat_thresh_high, brighten_sat_pwr, brighten_magnitude, new_aces)*/ \
      brighten_dark_sat(aces_d, .2/*ooo.175*//*.25*/, 0.0/*ooo.275*/, .85/*ooo1.0*/, .8/*ooo.9*/, .16/*jjj.13*//*hhh0.20*//*ooo.35*//*ooo.325*//*.28*//*.24*//*.225*/,  aces_e) /* brighten dark saturated colors */ \
/*#define boost_sat(aces, sat_pwr, sat_range, bright_pwr, bright_max, sat_boost, new_aces) */ \
      boost_sat(aces_e, 1.7, .6, BRIGHT_PWR_FOR_SAT_BOOST, BRIGHT_MAX_FOR_SAT_BOOST, SAT_BOOST_AMOUNT, aces_f) /* boost saturation for middle-bright */ \
      boost_sat(aces_f, 1.0/*ooo1.4*/, .7, .85/*ooo.85*/, .2, .2/*jjj.275*//*ooo.45*//*mmm.325*/, aces_full) /* boost saturation in shadows */ \
   }


#define BRIGHT_REDUCE 0.65 // rolloff factor, should be less than one, designed to reduce brights */
#define BRIGHT_ASYMP_HALF 2.5 // half-way point for asymptote rolloff */

#ifdef BRIGHT_HIGHLIGHTS_IN_TONE_CURVE

 #define roll_off_brights(out, in) /* null operation */ \
  out[0]=in[0]; out[1]=in[1]; out[2]=in[2];

#else /* not BRIGHT_HIGHLIGHTS_IN_TONE_CURVE */
 #define roll_off_brights(out, in) \
{ \
      float rgb_norm, norm_asymp, weight, in_pos[3]; \
\
      in_pos[0] = MAX(TINY, in[0]); \
      in_pos[1] = MAX(TINY, in[1]); \
      in_pos[2] = MAX(TINY, in[2]); \
\
      weighted_yellow_third_norm(in_pos, rgb_norm); \
\
      norm_asymp = rgb_norm / (rgb_norm + BRIGHT_ASYMP_HALF); /* asymptote with BRIGHT_ASYMP_HALF as half-way point */ \
      weight = (1.0 - norm_asymp) + norm_asymp * BRIGHT_REDUCE; \
\
      out[0] = in[0] * weight; \
      out[1] = in[1] * weight; \
      out[2] = in[2] * weight; \
}
#endif /* BRIGHT_HIGHLIGHTS_IN_TONE_CURVE or not */

/************************************************************************************************************/
 #ifdef BYPASS_LMT
  #define lmt_position_processing(aces_in, aces_adj) { \
    aces_adj[0] = aces_in[0];  aces_adj[1] = aces_in[1];  aces_adj[2] = aces_in[2]; \
  }
 #else /* not BYPASS_LMT */
  #ifdef SIMPLE_LMT
   #define lmt_position_processing(aces_in, aces_adj) \
     simple_processing(aces_in, aces_adj)

  #else /* NOT SIMPLE_LMT */
   #ifdef GAMMA_AND_MAT
    #define lmt_position_processing(aces_in, aces_adj) \
      gamma_plus_mat(aces_in, aces_adj)

   #else /* not GAMMA_AND_MAT */

    #ifdef FULL_LMT
#define OFF_DIAG_FULL 0.0/*jjj.01*//*jjj-.006*/ /* global saturation adjusted (postive is desat, negative is sat boost) */
     #define lmt_position_processing(aces_in, aces_adj) \
     { \
       full_processing(aces_in, aces_adj) \
       float off_diag = OFF_DIAG_FULL; \
       float diag = 1.0 - off_diag - off_diag; \
       float rr, gg, bb; \
       rr = aces_adj[0] * diag + MAX(0.0, aces_adj[1]) * off_diag \
                               + MAX(0.0, aces_adj[2]) * off_diag; \
       gg = aces_adj[1] * diag + MAX(0.0, aces_adj[0]) * off_diag \
                               + MAX(0.0, aces_adj[2]) * off_diag; \
       bb = aces_adj[2] * diag + MAX(0.0, aces_adj[0]) * off_diag \
                               + MAX(0.0, aces_adj[1]) * off_diag; \
       aces_adj[0] = rr; \
       aces_adj[1] = gg; \
       aces_adj[2] = bb; \
       float aces_rolloff[3]; \
       roll_off_brights(aces_rolloff, aces_adj) \
       aces_adj[0] = aces_rolloff[0]; \
       aces_adj[1] = aces_rolloff[1]; \
       aces_adj[2] = aces_rolloff[2]; \
     }
    #else /* not FULL_LMT */
     #ifdef DOUBLE_LMT
#define BLEND_PROPORTION 0.65 /* greater is more full, less is more moderate */
#define OFF_DIAG_DOUBLE -.0125 /* global saturation boost */
       #define lmt_position_processing(aces_in, aces_adj) \
       { float aces_adj_full[3]; \
         full_processing(aces_in, aces_adj_full) \
\
         moderate_processing(aces_in, aces_adj) \
\
         aces_adj[0] = BLEND_PROPORTION * aces_adj_full[0] + (1.0 - BLEND_PROPORTION) * aces_adj[0]; \
         aces_adj[1] = BLEND_PROPORTION * aces_adj_full[1] + (1.0 - BLEND_PROPORTION) * aces_adj[1]; \
         aces_adj[2] = BLEND_PROPORTION * aces_adj_full[2] + (1.0 - BLEND_PROPORTION) * aces_adj[2]; \
\
         float off_diag = OFF_DIAG_DOUBLE; \
         float diag = 1.0 - off_diag - off_diag; \
         float rr, gg, bb; \
         rr = aces_adj[0] * diag + MAX(0.0, aces_adj[1]) * off_diag \
                                 + MAX(0.0, aces_adj[2]) * off_diag; \
         gg = aces_adj[1] * diag + MAX(0.0, aces_adj[0]) * off_diag \
                                 + MAX(0.0, aces_adj[2]) * off_diag; \
         bb = aces_adj[2] * diag + MAX(0.0, aces_adj[0]) * off_diag \
                                 + MAX(0.0, aces_adj[1]) * off_diag; \
         aces_adj[0] = rr; /* note: these can be pushed negative if fully saturated near the ACES gamut edge */ \
         aces_adj[1] = gg; \
         aces_adj[2] = bb; \
       }
     #else /* not DOUBLE_LMT */
      #ifdef MODERATE_LMT
       #define lmt_position_processing(aces_in, aces_adj) \
         moderate_processing(aces_in, aces_adj)
      #else /* not MODERATE_LMT */

 error, no processing configuration selected

      #endif /* MODERATE_LMT or not */
     #endif /* DOUBLE_LMT or not */
    #endif /* FULL_LMT or not */
   #endif /* GAMMA_AND_MAT or not */
  #endif /* SIMPLE_LMT or not */
 #endif /* BYPASS_LMT or not */



/************************************************************************************************************/
#define RRT_GD10(aces, oces) \
  { \
\
float aces_adjusted[3]; \
      lmt_position_processing(aces, aces_adjusted) \
\
      float rgb_norm, rgb_norm_Post; \
\
      weighted_yellow_third_norm(aces_adjusted, rgb_norm); \
\
      rrt_shaper_fwd_gd10( rgb_norm, rgb_norm_Post) /* lookup using norm and spline-interpolate */ \
\
      if (rgb_norm > TINY_FOR_NORMS) { \
/* scale by the ratio of the norm to the looked-up-and-interpolated tone curve value of that norm */ \
/* This is RGB ratio preserving (thus chromaticity preserving).  This also preserves negative numbers in ratio proportion to positive numbers. */ \
        oces[0] = (aces_adjusted[0] * rgb_norm_Post) / rgb_norm; /* (rgb_norm > TINY_FOR_NORMS) prevents divide by zero */ \
        oces[1] = (aces_adjusted[1] * rgb_norm_Post) / rgb_norm; \
        oces[2] = (aces_adjusted[2] * rgb_norm_Post) / rgb_norm; \
      } else { /* Below TINY_FOR_NORMS, unity rrt_shaper_fwd and ratio restore allow values to be simply scaled by .02 slope. */ \
        /* This also passes (scaled) negative numbers. */ \
        /* note that the key ingredient to make this continuous is matched .02 slope in rrt_shaper_fwd for small values of rgb_norm (below -3.5 log DARK_START corresponding to TINY_FOR_NORMS) */ \
        oces[0] = .02 * aces_adjusted[0]; \
        oces[1] = .02 * aces_adjusted[1]; \
        oces[2] = .02 * aces_adjusted[2]; \
      } /* rgb_norm > TINY_FOR_NORMS or not */ \
      oces[0] = MAX(TINY, oces[0]); \
      oces[1] = MAX(TINY, oces[1]); \
      oces[2] = MAX(TINY, oces[2]); \
} /* RRT_GD10 */






