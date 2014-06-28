
/* sigma compare, author Gary Demos June 2014 */
/* no warranties expressed nor implied */
/* no representation is herein made as to usefulness nor suitability for any pupose */
/* caution, code may contain bugs and design flaws */
/* use at your own risk */

// to build:
/* 

g++ sigma_compare.cpp dpx_file_io.cpp -o sigma_compare_linux -lpthread \
 -I /usr/local/include/OpenEXR \
 /usr/local/lib/libHalf.a /usr/local/lib/libIlmImf.a /usr/local/lib/libIex.a /usr/local/lib/libIlmThread.a \
 /usr/local/lib/libz.a -m64 -O1

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include <OpenEXR/ImfRgbaFile.h> /* OpenExr */
#include <OpenEXR/ImfArray.h> /* OpenExr */

using namespace Imf; /* OpenExr */
using namespace Imath; /* OpenExr */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define NUMBER_OF_STOPS 32 // note: the 5 bit half-float exponent is limited to a range of 31 stops (value 31 = +-inf nan out-of-range code).  Note: value 0-exp denorm extends lower),  Use 32-bit floating point I/O (e.g. dpx32) for ranges greater than 31 (+denorm) stops

// Prototypes:
void dpx_read (char *inname, float **pixels_read, short *width, short *height, short cineon, short half_flag);



/***********************************************************************************************************/
int
main(int argc, char **argv)
{

short x, y, c, j, k;
float *pixels = NULL, *pixels2 = NULL, *pixels3 = NULL;
short *bin_value = NULL;
char infile1[300], infile2[300];

int iframe=0, first_frame=0, last_frame=0;

short h_reso, v_reso;
float col[3];
float flt_log2;
double square_sum_col[3][NUMBER_OF_STOPS];
double average_col[3][NUMBER_OF_STOPS];
float sigma_col[3][NUMBER_OF_STOPS];
double furthest_outlier_avg[3][NUMBER_OF_STOPS];

double count_col[3][NUMBER_OF_STOPS]; /* must be double precision 64-bit floating point to continue to accumuate counts beyond the 23-bit lsb (could also use long-long 64-bit int, but 32-bit int is insufficient) */
#define SIGMA_MULTIPLES 7 /*2,3,4,6,8,12,16*/
double count_col_multiple_sigma[SIGMA_MULTIPLES][3][NUMBER_OF_STOPS];
double count_neg_col[3];
float col_dif[3];
double average_neg_col[3];
double square_sum_neg_col[3];
float sigma_neg_col[3];
short dpx_in_file1=0, dpx_in_file2=0;
float amplitude_factor;

    Array2D<Rgba> half_float_pixels; /* unused if both input files are dpx */

if (argc <= 2) {
  printf(" usage: %s original_files, test_comparison_files, first_frame, last_frame, amplitude_factor\n", argv[0]);
  printf(" note: use %% format for frame numbers, such as %%07d, within the file names\n");
  exit(1);
}

    if (argc > 4) {
      first_frame = atoi(argv[3]);
      last_frame = atoi(argv[4]);
      printf(" first_frame = %d last_frame = %d\n", first_frame, last_frame);

    } else { /* argc <= 4 */
      printf(" no frame range, one frame will be processed\n");
    } /* argc > 4 or not */

    if (argc > 5) {
      amplitude_factor = atof(argv[5]);
      printf(" setting amplitude_factor to %f\n", amplitude_factor);
 /* can scale nits to unity relationship here
 e.g. if 1.0=100nits, use 100.0 here for the amplitude_factor in order to see the printout relative to nits (sigma and average value), does not affect self_relative
 e.g. also if looking at nits, but wish to see in terms of 100nits = 1.0, then use .01 for amplitude_factor to range 100nits down to 1.0 (e.g. a nominal white at 1.0) */
    } else { /* argc <= 5 */
      amplitude_factor = 1.0;
      printf(" using default amplitude_factor = %f\n", amplitude_factor);
    } /* argc > 5 or not */


    for (c=0; c<3; c++) {
      for (j=0; j<NUMBER_OF_STOPS; j++) {
        square_sum_col[c][j] = 0.0;
        average_col[c][j] = 0.0;
        count_col[c][j] = 0;
        furthest_outlier_avg[c][j] = 0;
        for(k=0; k < SIGMA_MULTIPLES; k++) {
          count_col_multiple_sigma[k][c][j] = 0;
        } /* k */
      } /* j */
      square_sum_neg_col[c] = 0.0;
      count_neg_col[c] = 0;
      average_neg_col[c] = 0;
    } /* c */


  for (iframe = first_frame; iframe <= last_frame; iframe++) {
    sprintf(infile1, argv[1], iframe);
    sprintf(infile2, argv[2], iframe);

    short num_chars = strlen(infile1); /* length of in_file_name string */

    if ((!strcmp(&infile1[num_chars-1],  "x"))  ||(!strcmp(&infile1[num_chars-1],    "X")) || /* DPX file ending in ".dpx" (not recommended, 10bit linear has insufficient precision near black) */
        (!strcmp(&infile1[num_chars-3],  "x32"))||(!strcmp(&infile1[num_chars-3],  "X32")) || /* dpx32 or DPX32 */
        (!strcmp(&infile1[num_chars-3],  "x16"))||(!strcmp(&infile1[num_chars-3],  "X16")) || /* dpx16 or DPX16 (not recommended) */
        (!strcmp(&infile1[num_chars-4], "xhlf"))||(!strcmp(&infile1[num_chars-4], "XHLF"))) { /* dpxhlf or DPXHLF, 16-bit half-float (non-standard, improperly defined in dpx 2.0 standard) */
      dpx_in_file1 = 1;
    }
    num_chars = strlen(infile2); /* length of in_file_name string */

    if ((!strcmp(&infile2[num_chars-1],  "x"))  ||(!strcmp(&infile2[num_chars-1],    "X")) || /* DPX file ending in ".dpx" (not recommended, 10bit linear has insufficient precision near black) */
        (!strcmp(&infile2[num_chars-3],  "x32"))||(!strcmp(&infile2[num_chars-3],  "X32")) || /* dpx32 or DPX32 */
        (!strcmp(&infile2[num_chars-3],  "x16"))||(!strcmp(&infile2[num_chars-3],  "X16")) || /* dpx16 or DPX16 (not recommended) */
        (!strcmp(&infile2[num_chars-4], "xhlf"))||(!strcmp(&infile2[num_chars-4], "XHLF"))) { /* dpxhlf or DPXHLF, 16-bit half-float (non-standard, improperly defined in dpx 2.0 standard) */
      dpx_in_file2 = 1;
    }


    if (dpx_in_file2 == 1) {

      dpx_read (infile2, &pixels, &h_reso, &v_reso, 0/*not cineon*/, 0/*no half_flag*/);

      if (amplitude_factor != 1.0) {

        for (c=0; c<3; c++) {
          for(y=0; y< v_reso; y++) {
            for(x=0; x< h_reso; x++) {
              pixels[(c*v_reso + y) * h_reso + x] = amplitude_factor * pixels[(c*v_reso + y) * h_reso + x];
            } /* x */
          } /* y */
        } /* c */
      } /* amplitude_factor != 1.0 */

    } else { /* exr */

{
    RgbaInputFile file (infile2, 1);

    Box2i dw = file.dataWindow();

    h_reso = dw.max.x - dw.min.x + 1;
    v_reso = dw.max.y - dw.min.y + 1;
    half_float_pixels.resizeErase (v_reso, h_reso);

    file.setFrameBuffer (&half_float_pixels[0][0] - dw.min.x - dw.min.y * h_reso, 1, h_reso);
    file.readPixels (dw.min.y, dw.max.y);

    printf(" reading target comparison file %s having h_reso = %d v_reso = %d\n", infile2, h_reso, v_reso);

    if (pixels == NULL) { pixels = (float *) malloc(h_reso * v_reso * 12); /* 4-bytes/float * 3-colors */ }

    for(y=0; y< v_reso; y++) {
      for(x=0; x< h_reso; x++) {
        pixels[(0*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].r;
        pixels[(1*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].g;
        pixels[(2*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].b;
      } /* x */
    } /* y */
}
    } /* dpx vs exr */
/**********************************************************************************************/

    if (dpx_in_file1 == 1) {
short h_reso_b, v_reso_b;

      dpx_read (infile1, &pixels3, &h_reso_b, &v_reso_b, 0/*not cineon*/, 0/*no half_flag*/);

      if ((h_reso_b != h_reso) || (v_reso_b != v_reso)) {
        printf(" error, file resolutions do not match, infile1 h_reso = %d v_reso = %d, infile2 h_reso = %d v_reso = %d, aborting\n",
          h_reso, v_reso, h_reso_b, v_reso_b);
        exit(1);
      }

      if (amplitude_factor != 1.0) {
        for (c=0; c<3; c++) {
          for(y=0; y< v_reso; y++) {
            for(x=0; x< h_reso; x++) {
              pixels3[(c*v_reso + y) * h_reso + x] = amplitude_factor * pixels3[(c*v_reso + y) * h_reso + x];
            } /* x */
          } /* y */
        } /* c */
      } /* amplitude_factor != 1.0 */

    } else { /* exr */

{
    RgbaInputFile file (infile1, 1);

    Box2i dw = file.dataWindow();

short h_reso_b, v_reso_b;
    h_reso_b = dw.max.x - dw.min.x + 1;
    v_reso_b = dw.max.y - dw.min.y + 1;

    if ((h_reso_b != h_reso) || (v_reso_b != v_reso)) {
      printf(" error, file resolutions do not match, infile1 h_reso = %d v_reso = %d, infile2 h_reso = %d v_reso = %d, aborting\n",
        h_reso, v_reso, h_reso_b, v_reso_b);
      exit(1);
    }

    half_float_pixels.resizeErase (v_reso, h_reso);

    file.setFrameBuffer (&half_float_pixels[0][0] - dw.min.x - dw.min.y * h_reso, 1, h_reso);
    file.readPixels (dw.min.y, dw.max.y);

    printf(" reading original file %s having h_reso = %d v_reso = %d\n", infile1, h_reso, v_reso);
}
  } /* dpx vs exr file1 */

    if (pixels2 == NULL) { pixels2 = (float *) malloc(h_reso * v_reso * 12); /* 4-bytes/float * 3-colors */ }
    if (bin_value == NULL) { bin_value = (short *) malloc(h_reso * v_reso * 6); /* 2-bytes/short * 3-colors */ }


    for (c=0; c<3; c++) {
      for(y=0; y< v_reso; y++) {
        for(x=0; x< h_reso; x++) {
          if (dpx_in_file1 == 1) {
            col[c] = pixels3[(c*v_reso + y) * h_reso + x];
          } else { /* exr */
            if (c==0) { col[c] = amplitude_factor * half_float_pixels[y][x].r; }
            if (c==1) { col[c] = amplitude_factor * half_float_pixels[y][x].g; }
            if (c==2) { col[c] = amplitude_factor * half_float_pixels[y][x].b; }
          } /* dpx vs exr */
/* store difference from original in pixels2 */
          col_dif[c] = col[c] - pixels[(c*v_reso + y) * h_reso + x]; 
          pixels2[(c*v_reso + y) * h_reso + x] = col_dif[c];

          if (col[c] < 0.0) {
            square_sum_neg_col[c] = square_sum_neg_col[c] + col_dif[c] * col_dif[c];
            average_neg_col[c] = average_neg_col[c] + col[c];
            count_neg_col[c] = count_neg_col[c] + 1;
            bin_value[(c*v_reso + y) * h_reso + x] = -1;
          } else { /* col[c] >= 0 */
            flt_log2 = log2f(MAX(1e-16, col[c]));
            j = flt_log2; /* float to int */
            if (flt_log2 > 0.0) j = j + 1; /* must round one side of zero for continuity of integer j */
            j = MIN(NUMBER_OF_STOPS - 1, MAX(0, j + NUMBER_OF_STOPS / 2)); /* split NUMBER_OF_STOPS half above 1.0 and half below 1.0 */
            square_sum_col[c][j] = square_sum_col[c][j] + col_dif[c] * col_dif[c];
            average_col[c][j] = average_col[c][j] + col[c];
            count_col[c][j] = count_col[c][j] + 1;
            bin_value[(c*v_reso + y) * h_reso + x] = j;
          } /* col[c] < 0 or not */
        } /* x */
      } /* y */
    } /* c */

    if (dpx_in_file1 == 1) {
      free(pixels3);
    } /* dpx */
    if (dpx_in_file2 == 1) {
      free(pixels);
    } /* dpx */

  } /* iframe */

  printf("\n amplitude_factor = %f\n", amplitude_factor);

    for (c=0; c<3; c++) {
printf("\n");
      if (count_neg_col[c] > 0) {
        square_sum_neg_col[c] = square_sum_neg_col[c] / count_neg_col[c];
        sigma_neg_col[c] = sqrtf(square_sum_neg_col[c]);
        average_neg_col[c] = average_neg_col[c] / count_neg_col[c];
        if (c==0) { printf(" sigma_neg_red = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
        if (c==1) { printf(" sigma_neg_grn = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
        if (c==2) { printf(" sigma_neg_blu = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
      } /* count_neg_col[c] > 0 */

      for (j=0; j<NUMBER_OF_STOPS; j++) {
        if (count_col[c][j] > 0) {
          square_sum_col[c][j] = square_sum_col[c][j] / count_col[c][j];
          sigma_col[c][j] = sqrtf(square_sum_col[c][j]);
          average_col[c][j] = average_col[c][j]  / count_col[c][j];
//          if (j == 0) { average_col[c][j] = MAX(1e-16, average_col[c][j]); /* prevent divide by zero */ }
          if (c==0) { printf(" sigma_red[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==1) { printf(" sigma_grn[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==2) { printf(" sigma_blu[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
        } /* count_col[c][j] > 0 */
      } /* j */
    } /* c */

printf("\n beginning second pass over frames (multiples of sigma, now that sigma for all frames has been determined)\n\n");

    for (iframe = first_frame; iframe <= last_frame; iframe++) {
      sprintf(infile1, argv[1], iframe);
      sprintf(infile2, argv[2], iframe);


    if (dpx_in_file2 == 1) {

      dpx_read (infile2, &pixels, &h_reso, &v_reso, 0/*not cineon*/, 0/*no half_flag*/);

      if (amplitude_factor != 1.0) {

        for (c=0; c<3; c++) {
          for(y=0; y< v_reso; y++) {
            for(x=0; x< h_reso; x++) {
              pixels[(c*v_reso + y) * h_reso + x] = amplitude_factor * pixels[(c*v_reso + y) * h_reso + x];
            } /* x */
          } /* y */
        } /* c */
      } /* amplitude_factor != 1.0 */

    } else { /* exr */

{
    RgbaInputFile file (infile2, 1);

    Box2i dw = file.dataWindow();

    half_float_pixels.resizeErase (v_reso, h_reso);

    file.setFrameBuffer (&half_float_pixels[0][0] - dw.min.x - dw.min.y * h_reso, 1, h_reso);
    file.readPixels (dw.min.y, dw.max.y);

    printf(" reading target comparison file %s having h_reso = %d v_reso = %d\n", infile2, h_reso, v_reso);

    for(y=0; y< v_reso; y++) {
      for(x=0; x< h_reso; x++) {
        pixels[(0*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].r;
        pixels[(1*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].g;
        pixels[(2*v_reso + y) * h_reso + x] = amplitude_factor * half_float_pixels[y][x].b;
      } /* x */
    } /* y */
}
    } /* dpx vs exr */
/**********************************************************************************************/

    if (dpx_in_file1 == 1) {
short h_reso_b, v_reso_b;

      dpx_read (infile1, &pixels3, &h_reso_b, &v_reso_b, 0/*not cineon*/, 0/*no half_flag*/);

      if (amplitude_factor != 1.0) {
        for (c=0; c<3; c++) {
          for(y=0; y< v_reso; y++) {
            for(x=0; x< h_reso; x++) {
              pixels3[(c*v_reso + y) * h_reso + x] = amplitude_factor * pixels3[(c*v_reso + y) * h_reso + x];
            } /* x */
          } /* y */
        } /* c */
      } /* amplitude_factor != 1.0 */

    } else { /* exr */

{
    RgbaInputFile file (infile1, 1);

    Box2i dw = file.dataWindow();

    half_float_pixels.resizeErase (v_reso, h_reso);

    file.setFrameBuffer (&half_float_pixels[0][0] - dw.min.x - dw.min.y * h_reso, 1, h_reso);
    file.readPixels (dw.min.y, dw.max.y);

    printf(" reading original file %s having h_reso = %d v_reso = %d\n", infile1, h_reso, v_reso);
}
  } /* dpx vs exr file1 */


/* recompute pixels2 and bin_value */
    for (c=0; c<3; c++) {
      for(y=0; y< v_reso; y++) {
        for(x=0; x< h_reso; x++) {
          if (dpx_in_file1 == 1) {
            col[c] = pixels3[(c*v_reso + y) * h_reso + x];
          } else { /* exr */
            if (c==0) { col[c] = amplitude_factor * half_float_pixels[y][x].r; }
            if (c==1) { col[c] = amplitude_factor * half_float_pixels[y][x].g; }
            if (c==2) { col[c] = amplitude_factor * half_float_pixels[y][x].b; }
          } /* dpx vs exr */
/* store difference from original in pixels2 */
          col_dif[c] = col[c] - pixels[(c*v_reso + y) * h_reso + x]; 
          pixels2[(c*v_reso + y) * h_reso + x] = col_dif[c];

          if (col[c] < 0.0) {
            bin_value[(c*v_reso + y) * h_reso + x] = -1;
          } else { /* col[c] >= 0 */
            flt_log2 = log2f(MAX(1e-16, col[c]));
            j = flt_log2; /* float to int */
            if (flt_log2 > 0.0) j = j + 1; /* must round one side of zero for continuity of integer j */
            j = MIN(NUMBER_OF_STOPS - 1, MAX(0, j + NUMBER_OF_STOPS / 2)); /* split NUMBER_OF_STOPS half above 1.0 and half below 1.0 */
            bin_value[(c*v_reso + y) * h_reso + x] = j;
          } /* col[c] < 0 or not */
        } /* x */
      } /* y */
    } /* c */


    for (c=0; c<3; c++) {
      for(y=0; y< v_reso; y++) {
        for(x=0; x< h_reso; x++) {

          if (bin_value[(c*v_reso + y) * h_reso + x] >= 0) { /* dont count negs */
            j = bin_value[(c*v_reso + y) * h_reso + x];
            if (count_col[c][j] > 0) {
              col_dif[c] = pixels2[(c*v_reso + y) * h_reso + x];

              for(k=0; k < SIGMA_MULTIPLES; k++) {
short multiple;
                multiple = (1 << ((k>>1) + 1));
                multiple = multiple + (k & 1) * (multiple >> 1); /* turns k=0,1,2,3,4,5,6,7 into multiple=2,3,4,6,8,12,16 */
                if (fabsf(col_dif[c]) > (multiple * sigma_col[c][j])) { count_col_multiple_sigma[k][c][j] = count_col_multiple_sigma[k][c][j] + 1;

//if((k== (SIGMA_MULTIPLES-1))&&(j==16)){printf(" c=%d x=%d y=%d col_dif[%d]=%e\n", c, x, y, c, col_dif[c]);}
                furthest_outlier_avg[c][j] = furthest_outlier_avg[c][j] + fabsf(col_dif[c]);
 }
              } /* k */
            } /* count_col[c][j] > 0 */
          } /* bin_value[(c*v_reso + y) * h_reso + x] >= 0 */

        } /* x */
      } /* y */
    } /* c */

    if (dpx_in_file1 == 1) {
      free(pixels3);
    } /* dpx */
    if (dpx_in_file2 == 1) {
      free(pixels);
    } /* dpx */

   } /* iframe */


   printf("\n amplitude_factor = %f\n", amplitude_factor);

    for (c=0; c<3; c++) {
printf("\n");
      if (count_neg_col[c] > 0) {
        if (c==0) { printf(" sigma_neg_red = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
        if (c==1) { printf(" sigma_neg_grn = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
        if (c==2) { printf(" sigma_neg_blu = %e self_relatve = %f (%f%%) at average value = %e for %.0f pixels\n",
          sigma_neg_col[c], sigma_neg_col[c] / average_neg_col[c], 100.0 *sigma_neg_col[c] / average_neg_col[c], average_neg_col[c], count_neg_col[c]); }
      } /* count_neg_col[c] > 0 */


      for (j=0; j<NUMBER_OF_STOPS; j++) {
        if (count_col[c][j] > 0) {
          if (c==0) { printf(" sigma_red[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==1) { printf(" sigma_grn[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==2) { printf(" sigma_blu[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
        } /* count_col[c][j] > 0 */
      } /* j */
    } /* c */


short multiple;
    for (c=0; c<3; c++) {
printf("\n");
      for (j=0; j<NUMBER_OF_STOPS; j++) {
        if ((count_col[c][j] > 0) && (count_col_multiple_sigma[0][c][j] > 0)) {
          if (c==0) { printf("  sigma_red[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==1) { printf("  sigma_grn[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
          if (c==2) { printf("  sigma_blu[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels\n",
            j, sigma_col[c][j], sigma_col[c][j] / average_col[c][j], 100.0 * sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col[c][j]); }
        } /* count_col[c][j] > 0 */
        for(k=0; k < (SIGMA_MULTIPLES-1); k++) {
          if ((count_col_multiple_sigma[k][c][j] > 0) && (count_col_multiple_sigma[k][c][j] > count_col_multiple_sigma[k+1][c][j])) {
            multiple = (1 << ((k>>1) + 1));
            multiple = multiple + (k & 1) * (multiple >> 1); /* turns k=0,1,2,3,4,5,6,7 into multiple=2,3,4,6,8,12,16 */
            if (c==0) { printf(" %dsigma_red[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, multiple*sigma_col[c][j], multiple*sigma_col[c][j] / average_col[c][j], 100.0 * multiple*sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
            if (c==1) { printf(" %dsigma_grn[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, multiple*sigma_col[c][j], multiple*sigma_col[c][j] / average_col[c][j], 100.0 * multiple*sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
            if (c==2) { printf(" %dsigma_blu[%d] = %e self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, multiple*sigma_col[c][j], multiple*sigma_col[c][j] / average_col[c][j], 100.0 * multiple*sigma_col[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
          } /* count_col_multiple_sigma[k][c][j] > 0 */
        } /* k */
        k = SIGMA_MULTIPLES-1;
        if (count_col_multiple_sigma[k][c][j] > 0) {
            furthest_outlier_avg[c][j] = furthest_outlier_avg[c][j] / count_col_multiple_sigma[k][c][j]; /* compute average from sum */
            multiple = (1 << ((k>>1) + 1));
            multiple = multiple + (k & 1) * (multiple >> 1); /* turns k=0,1,2,3,4,5,6,7 into multiple=2,3,4,6,8,12,16 */
            if (c==0) { printf(" beyond %dsigma_red[%d] furthest_outlier_avg dif = %e which is %f sigma (sigma being %e), self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, furthest_outlier_avg[c][j], furthest_outlier_avg[c][j] / sigma_col[c][j], sigma_col[c][j], furthest_outlier_avg[c][j] / average_col[c][j], 100.0 * furthest_outlier_avg[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
            if (c==1) { printf(" beyond %dsigma_grn[%d] furthest_outlier_avg dif = %e which is %f sigma (sigma being %e), self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, furthest_outlier_avg[c][j], furthest_outlier_avg[c][j] / sigma_col[c][j], sigma_col[c][j], furthest_outlier_avg[c][j] / average_col[c][j], 100.0 * furthest_outlier_avg[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
            if (c==2) { printf(" beyond %dsigma_blu[%d] furthest_outlier_avg dif = %e which is %f sigma (sigma being %e), self_relative = %f (%f%%) at average value = %e for %.0f pixels which is %f %% of the %.0f pixels within a stop of this value\n",
              multiple, j, furthest_outlier_avg[c][j], furthest_outlier_avg[c][j] / sigma_col[c][j], sigma_col[c][j], furthest_outlier_avg[c][j] / average_col[c][j], 100.0 * furthest_outlier_avg[c][j] / average_col[c][j], average_col[c][j], count_col_multiple_sigma[k][c][j], (count_col_multiple_sigma[k][c][j] * 100.0)/(1.0 * count_col[c][j]), count_col[c][j]); }
        } /* (count_col_multiple_sigma[k][c][j] > 0) */
      } /* j */
    } /* c */


} /* main */


