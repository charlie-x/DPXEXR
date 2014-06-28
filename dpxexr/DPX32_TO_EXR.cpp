

/* ///////////////////////////////////////////////////////////////////////////
//
//
// Copyright (c) 2014 Gary Demos
// Author: Gary Demos
//
//////////////////////////////////////////////////////////////////////////////

/* no warranties expressed nor implied */
/* no representation is herein made as to usefulness nor suitability for any pupose */
/* caution, code may contain bugs and design flaws */
/* use at your own risk */

//----------------------------------------------------------------------------
//
// DPX32_TO_EXR.cpp
// Intended for function on m64 (64-bit OS) MacOSX Snow Leopard (m32 should work on MacOSX as well) and x86_64 (64-bit OS) Linux 
//


// usage: DPX32_TO_EXR_xxx infiles, outfiles, first_frame, last_frame, first_frame_out


// to build:
// linux:
/*
g++ DPX32_TO_EXR.cpp dpx_file_io.cpp -o DPX32_TO_EXR_LINUX -lpthread \
 -I /usr/local/include/OpenEXR \
 /usr/local/lib/libHalf.a /usr/local/lib/libIlmImf.a /usr/local/lib/libIex.a /usr/local/lib/libIlmThread.a \
 /usr/local/lib/libz.a -m64

*/
// MacOSX:
/*
g++ DPX32_TO_EXR.cpp dpx_file_io.cpp -o DPX32_TO_EXR_MACOSX -lpthread \
 -I /usr/local/include/OpenEXR \
 /usr/local/lib/libHalf.a /usr/local/lib/libIlmImf.a /usr/local/lib/libIex.a /usr/local/lib/libIlmThread.a \
 /usr/local/lib/libz.a -m64

*/

// note: -O3 will make the code run faster, but some compilers generate bad executables with -O3, -O2, or even -O1.  Most compilers default the above to -O1

// OpenExr must be installed for libIlmImf library, see OpenExr developers if there are problems installing OpenExr
// zlib (libz) must also be installed for OpenEXR

// OpenExr headers included here (via ./OpenExr_includes/ subdirectory) correspond to OpenEXR-1.4.0 (old, but stable).  These headers should function on both MacOSX and Linux.
// These OpenExr headers have been modified to eliminate some includes within the headers, so that the headers can be kept self-contained within the OpenExr_includes subdirectory at this level
//
// See OpenExr.com website for OpenExr license information
//////////////////////////////////////////////////////////////////  */

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

#define INT_SW(A) (((A >> 24) & 255) | (((A >> 16) & 255) << 8) | (((A >> 8) & 255) << 16) | ((A & 255) << 24));
#define SHORT_SW(A) (((A >> 8) & 255) | ((A & 255) << 8));


#define INTEL_LE

// Prototypes:
void dpx_read (char *inname, float **pixels_read, short *width, short *height, short cineon, short half_flag);



  short h_reso, v_reso;

  Array2D<Rgba> half_float_pixels;

  float *pixels = NULL;

/***********************************************************************************************************/
int
main(int argc, char **argv)
{
short x,y,c,i;
char outfile[300], infile[300];
short num_chars;
int first, last, frame, first_out;
float s, t, tmp;
float red, grn, blu;
short argnm=0;



if (argc <= 4) {
  printf(" usage: %s infiles, outfiles, first_frame, last_frame, first_frame_out\n", argv[0]);
  exit(1);
}

 first = atoi(argv[3]);
 last  = atoi(argv[4]);
 if(argc > 5) {
   first_out = atoi(argv[5]);
 } else {
   first_out=first;
 }

 printf(" processing frames %d to %d to frames %d to %d\n", first, last, first_out, last + first_out - first);


/*******************************************************************************************************************************************************/
/* frame loop: */
   for (frame=first; frame <= last; frame++) {

     sprintf(infile, argv[1], frame);
     num_chars = strlen(infile); /* length of infile string */
     if ((!strcmp(&infile[num_chars-1], "x"))  ||(!strcmp(&infile[num_chars-1], "X")) || /* DPX floating point file ending in ".dpx" */
         (!strcmp(&infile[num_chars-3], "x32"))||(!strcmp(&infile[num_chars-3], "X32"))) { /* DPX32 */
       printf(" processing input file %s\n", infile);
     } else { /* not exr */
      printf(" unknown filetype for reading, since extension doesn't end in x or x32, only dpx float reading supported, infile = %s, aborting\n", outfile);
      exit(1);
     } /* exr or not */

     sprintf(outfile, argv[2], frame + first_out - first);
     num_chars = strlen(outfile); /* length of outfile string */
       if ((!strcmp(&outfile[num_chars-1], "r"))||(!strcmp(&outfile[num_chars-1], "R"))) { /* EXR file ending in ".exr" */
        printf(" processing output file %s\n", outfile);
       } else { /* not exr */
        printf(" unknown filetype for writing, since extension doesn't end in r, only exr writing supported, outfile = %s, aborting\n", outfile);
        exit(1);
       } /* exr or not */

    dpx_read (infile, &pixels, &h_reso, &v_reso, 0/*not cineon*/, 0/*no half_flag*/);

    printf(" begin writing %s with resolution y=%d x=%d\n", outfile, h_reso, v_reso);

    Array2D<Rgba> pixels_out_exr(v_reso, h_reso);

    for(c=0; c<3; c++) {
      for (y=0; y<v_reso; y++) {
        for (x=0; x<h_reso; x++) {
          tmp = pixels[(c * v_reso + y) * h_reso + x];
          if (c==0) pixels_out_exr[y][x].r = tmp;
          if (c==1) pixels_out_exr[y][x].g = tmp;
          if (c==2) pixels_out_exr[y][x].b = tmp;
        } /* x */
      } /* y */
    } /* c */

    free(pixels); pixels = NULL;

    RgbaOutputFile file_out (outfile, h_reso, v_reso, WRITE_RGB, 1 /* one thread for writing */ );

    file_out.setFrameBuffer (&pixels_out_exr[0][0], 1, h_reso);
    file_out.writePixels (v_reso);

    printf(" finished writing %s at x_reso = %d y_reso = %d\n", outfile, h_reso, v_reso);

 } /* frame loop */

} /* main */

