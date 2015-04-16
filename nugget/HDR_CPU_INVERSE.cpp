
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


/* ///////////////////////////////////////////////////////////////////////////
//
//
// Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, AMPAS
// Author: Gary Demos 2010,2011,2012,2013,2014, 2015
//
// Algorithms and code herein also are copyright Gary Demos, and some portions may retain a private use and copyright
//
// Further developed in collaboration with the ASC Technology Committee, beginning March 2015
//
//////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------
//
// HDR_CPU.cpp
// CPU example implementation (SDK) for HDR System
//
// This code inverts the process of going from ACES to display-referred dpx files.  The intent is that graded dpx files can be converted to pre-nugget half-floats for repurposing using forward processing


//
// to build:
// linux (tested on Scientific Linux 6.6):
g++ HDR_CPU_INVERSE.cpp -o HDR_CPU_INVERSE_LINUX -lpthread \
 -I /usr/local/include/OpenEXR \
 /usr/local/lib/libHalf.a /usr/local/lib/libIlmImf.a /usr/local/lib/libIex.a /usr/local/lib/libIlmThread.a \
 /usr/local/lib/libz.a -O3 -m64

// note: -O3 will make the code run faster, but some compilers generate bad executables with -O3, -O2, or even -O1.  Most compilers default to -O1

// OpenExr must be installed for libIlmImf library, see OpenExr developers if there are problems installing OpenExr
// zlib (libz) must also be installed for OpenExr

// No attempt is made in this code to read nor write OpenExr header information.
// The use of exr output files to contain pixels-ready-for-display does not indicate the ODT type nor ODT information in the exr header.

//
// Note: dpx read/write code provided by Gary Demos, retaining Gary Demos copyright 2010, 2011, 2012, 2013, 2014, 2015, no warranty implied, use or modify at your own risk

// No attempt is made in this code to write meaningful dpx header information (other than a generic header).
// Programs which rely on dpx header information may thus mis-interpret the generic dpx header which is written by this code.

// Note: the HDR System is intended for use with 32-bit float or OpenExr/GPU-style 16-bit half-float data in ACES RGB (having virtual R,G,B chromaticities)
// OpenExr headers included here (via ./OpenExr_includes/ subdirectory) correspond to OpenEXR-2.1.0 (./configure --disable-ilmbasetest)
// and ilmbase-2.1.0 and zlib-1.2.5.
// These OpenExr headers have been modified to eliminate some includes within the headers, so that the headers can be kept self-contained within the OpenExr_includes subdirectory at this level
//
// See OpenExr.com website for OpenExr license information
//
//----------------------------------------------------------------------------- */

// This implementation is an example, and does not contain the full range of intended capabilities nor intended alternate useful whitepoint processing possibilities.

// Gary Demos wishes to acknowledge the detailed help, clarifications, and algorithmic ideas provided by 
// M. Uchida, Iwaki, Ray Feeney, Jim Houston, Lars Borg, David Ruhoff, Jon McElvain, Jack Holm, Alex Forsythe, and Scott Dyer in the development of this sample C-code CPU implementation.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#include <OpenEXR/ImfRgbaFile.h> /* OpenExr */
#include <OpenEXR/ImfArray.h> /* OpenExr */

using namespace Imf; /* OpenExr */
using namespace Imath; /* OpenExr */

#include "hdr_macros.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define INT_SW(A) (((A >> 24) & 255) | (((A >> 16) & 255) << 8) | (((A >> 8) & 255) << 16) | ((A & 255) << 24));
#define SHORT_SW(A) (((A >> 8) & 255) | ((A & 255) << 8));

#define INTEL_LE

// globals:
short odt_through_flag = 0; /* 0 = normal odt, 1 = through odt, 2 = through odt with 0.1 scaledown */


#define LOOKUP_DIMENSION 256

#define RANGE_TOP_A .1 /* span 0.0 up to this amount with random numbers */
#define RANGE_TOP_B .6 /* span 0.0 up to this amount with random numbers */
#define RANGE_TOP_C 2.0 /* span 0.0 up to this amount with random numbers */
#define RANGE_TOP_D 12.0 /* span 0.0 up to this amount with random numbers */

/* also walk along the ACES gamut boundaries */
#define RNG_RR 32.0
#define RNG_GR 4.0
#define RNG_BR 4.0

#define RNG_RG 4.0
#define RNG_GG 32.0
#define RNG_BG 4.0

#define RNG_RB 4.0
#define RNG_GB 4.0
#define RNG_BB 32.0
/* span a large range */
#define RNG_LG 96.0
#define RNG_LG 96.0
#define RNG_LG 96.0

#define RANDOM_BATCH_SIZE 128
#define NUMBER_OF_RANDOM_BATCHES 64//256 /* 64 usually works, but may have a few bright speckle areas, 256 is a better choice (although slower) */

/* note: it is recommended that only the nugget (no LMT) be used, and that RADIOMETRIC_ODTS be active, with HDR/MDR_RADIOMETRIC_PROPORTION set to 1.0 */
/* this is expected to work with both HDR and MDR display-type transforms */

/* NOTE: THIS CODE IS EXPERIMENTAL AND COMPLETELY UNSUPPORTED */

short *count;
float *lookup[3];


// DPX file readers and writers
/***************************************************************************************************************************************************/
/* this makes the header 8192 */
#define USER_AREA_PAD 6144

static char  defaultDevice[] = "Unknown_Device";
static char  defaultNumber[] = "$Date: 2014/4/17 12:00:00 $";

typedef unsigned long U32;
typedef unsigned short U16;
typedef unsigned char U8;
typedef float R32;

typedef struct dpxFileHdr {
    U32		magicNumber;	  /* core */
    U32		imageOffset;	  /* core */
    char	version[8];	    /* core */
    U32		totalFileSize;	/* core */
    U32		dittoKey;
    U32		genericLength;
    U32		industrySpecificLength;
    U32		userLength;
    char	imageFileName[100];
    char	creationDateAndTime[24]; /* YYYY:MM:DD:HH:MM:SS:LTZ */
    char	creatorName[100];
    char	projectName[200];
    char	copyright[200];
    U32		encryptionKey;
    U8		reserved[104];
} dpxFileHdr;

typedef struct element {
    U32		dataSigned;	    /* core */
    U32		dataMin;
    R32		quantityMin;
    U32		dataMax;
    R32		quantityMax;
    U8		descriptor;	    /* core */
    U8		transfer;	      /* core */
    U8		colorimetric;	  /* core */
    U8		bitsPerElement;	/* core */
    U16		packing;	      /* core */
    U16		encoding;	      /* core */
    U32		offset;		      /* core */
    U32		endLinePadding;
    U32		endImagePadding;
    char	description[32];
} element;

typedef struct dpxImgHdr {
    U16		imageOrientation;     /* core */
    U16		numberOfElements;     /* core */
    U32		pixelsPerLine;        /* core */
    U32		linesPerElement;      /* core */
    element	elements[8];
    U8		reserved[52];
} dpxImgHdr;

typedef struct dpxOrientationHdr {
    U32		xOffset;
    U32		yOffset;
    R32		xCenter;
    R32		yCenter;
    U32		xOriginalSize;
    U32		yOriginalSize;
    char	originalImageFileName[100];
    char	originalDateAndTime[24];
    char	inputDevice[32];
    char	inputDeviceSerialNumber[32];
    U16		borderValidity[4];
    U32		pixelAspectRatio[2];
    U8		reserved[28];
} dpxOrientationHdr;

typedef struct dpxFilmHdr {
    char	filmManufacturingIdCode[2];
    char	filmType[2];
    char	perfsOffset[2];
    char	prefix[6];
    char	count[4];
    char	format[32];
    U32		framePosition;
    U32		frames;
    U32		heldCount;
    R32		frameRate;
    R32		shutterAngle;
    char	keyFrame[32];
    char	slateInfo[100];
    U8		reserved[56];
} dpxFilmHdr;

typedef struct dpxTVhdr {
    U32		timeCode;
    U32		userBits;
    U8		interlace;
    U8		fieldNumber;
    U8		videoSignalStandard;
    U8		unused;
    R32		horizontalFrequency;
    R32		verticalFrequency;
    R32		temporalFrequency;
    R32		syncToFirstPixelMicroSeconds;
    R32		gamma;
    R32		blackLevel;
    R32		blackGain;
    R32		breakpoint;
    R32		whiteLevel;
    R32		integrationTimeSeconds;
    U8		reserved[76];
} dpxTVhdr;

typedef struct dpxHdr {
    dpxFileHdr		    file;
    dpxImgHdr		    img;
    dpxOrientationHdr	src;
    dpxFilmHdr		    film;
    dpxTVhdr	        tv;
} DpxHdr;

static DpxHdr dpxHdr;
static int    dpxHdrInitialized = 0, Swapped = 0;

#define DPX_MAGIC		0x53445058   /* "SDPX" in ascii */
#define DPX_HEADER_SIZE		2048

#define UNDEF_U8		0xff
#define UNDEF_U16		0xffff
#define UNDEF_U32		0xffffffff
static U32 undefR32 = 0xffffffff;
#define UNDEF_R32 (*((float *)&undefR32))


#define CIN_MAGIC 0x802A5FD7
#define CIN_MAGIC_LE 0xD75F2A80

/***********************************************************************************************************/
/***********************************************************************************************************/
void
dpx_read_short(char *inname, unsigned short **pixels_read, short *width, short *height, short cineon)
{

FILE *fp;
unsigned short *pixels, *pixels_fr, *pxls_copy, *pxls_copy_fr;
short x,y;
unsigned int *pixelbuf;
unsigned int *pp;
unsigned int tmp, hdr_offset;
short wide, tall, ww, hh;
int cinhdr[64];
char dpx_header[2048];
short swap_byte_order;
int yvalr, yvalg, yvalb, yval1, yval2, yval3, yv1, yv2, yv3;

  swap_byte_order = 0;

  if ((fp = fopen(inname, "rb")) == NULL) {
    if (cineon != 0) {
      printf(" Cannot open cineon input file %s, aborting\n", inname);
    } else { /* dpx */
      printf(" Cannot open dpx input file %s, aborting\n", inname);
	} /* cineon vs dpx */
    exit(1);
  }

  if (cineon != 0) { /* cineon file */
  
    fread(cinhdr, 1, 256 /* 64 ints */, fp);
	
    tmp = cinhdr[0];
    if (tmp == CIN_MAGIC ) {
      printf(" reading cineon file big endian header \n ");
    } else {
      if (tmp == CIN_MAGIC_LE ) {
        printf(" reading cineon file little endian header \n ");
        swap_byte_order = 1;
	  } else { /* magic no good */
	    printf(" cineon magic number = %x no good, aborting\n", tmp);
		exit(1);
	  } /* cin_magic_le or not */
    } /* cin_magic (big endian) or not */


    if (swap_byte_order == 1) {
      hdr_offset = INT_SW(cinhdr[1]);
      tmp  = INT_SW(cinhdr[50]);
      wide = tmp; /* uint to short */
      tmp  = INT_SW(cinhdr[51]);
      tall = tmp; /* uint to short */
    } else { /* big_endian */
      hdr_offset = cinhdr[1];
      wide  = cinhdr[50]; /* uint to short */
      tall  = cinhdr[51]; /* uint to short */
    } /* swap_byte_order or not */

  } else { /* dpx file */

  fread( dpx_header, 1, 2048, fp); 

  if (ferror(fp)) {
    printf(" error reading header for file %s in dpx_read_short, aborting\n", inname);
    exit(1);
  }

  if ( *((int *) &dpx_header[0]) == 0x53445058) { /* big-endian "SDPX" in ascii */
    printf(" reading dpx big-endian file header \n ");
    swap_byte_order = 0;
  } else { /* not big-endian header */
    if ( *((int *) &dpx_header[0]) == 0x58504453) { /* little-endian "SDPX" in ascii */
      printf(" reading dpx little-endian file header \n ");
      swap_byte_order = 1;
    } else {
	  printf(" bad magic number = %x in dpx header read, aborting\n", *((int *) &dpx_header[0]));
      exit(1);
    } /* little-endian header or not */
  } /* big-endian header or not */

  if (swap_byte_order == 1) {
    tmp  = INT_SW(*((int *) &dpx_header[772]));
    wide = tmp; /* uint to short */
    tmp  = INT_SW(*((int *) &dpx_header[776]));
    tall = tmp; /* uint to short */
  } else { /* big_endian */
    wide = *((short *) &dpx_header[772]); /* uint to short */
    tall = *((short *) &dpx_header[776]); /* uint to short */
  } /* swap_byte_order or not */

  if (swap_byte_order == 1) {
    hdr_offset = INT_SW(*((unsigned int *) &dpx_header[4]));
  } else { /* big_endian */
    hdr_offset = *((unsigned int *) &dpx_header[4]);
  } /* swap_byte_order or not */

  if (dpx_header[803] == 10) {
#ifdef VERBOSE_DECODE
    printf(" dpx packing is normal 10-bit \n ");
#endif /* VERBOSE_DECODE */
  } else { /* not 10bit */
    printf(" dpx file is not 10-bit, but is rather %d bits, which are not support for dpx_read_short, aborting\n", dpx_header[803]);
    exit(1);
  } /* 10bit or not */
 } /* cineon vs. dpx header */


  printf(" width = %d, height = %d \n", wide, tall);
  *width = wide;
  *height = tall;

  fseek( fp, hdr_offset, SEEK_SET);
 
  if (ferror(fp)) {
    printf(" error seeking past header for file %s in dpx_read_short, aborting\n", inname);
    exit(1);
  }

  printf(" allocating memory for dpx pixels \n");

  pixelbuf = (unsigned int *) malloc(wide * tall * 4); /* 4 bytes per pixel in 10-bit dpx packing */
  pixels = (unsigned short *) malloc(wide * tall * 6); /* 2-bytes/short * 3-colors */
  *pixels_read = &pixels[0];

printf(" reading dpx file pixels \n ");

 for( y = 0; y < tall; y++) {

  pp = &pixelbuf[0] + y * wide;

  fread( pp, 1, wide * sizeof(unsigned int), fp);

  if (ferror(fp)) {
    printf(" error reading data for scanline %d for file %s in dpx_read_short, aborting\n", y, inname);
    exit(1);
  }

 } /* y */

 if(fclose(fp)) {
   printf(" error closing file %s in dpx_read_short, aborting\n", inname);
   exit(1);
 }

 for( y = 0; y < tall; y++) {
  yvalr = (0*tall + y)*wide;
  yvalg = (1*tall + y)*wide;
  yvalb = (2*tall + y)*wide;

  for( x = 0; x < wide ; x++) {


    if (swap_byte_order == 1) {
      tmp=INT_SW(pixelbuf[y * wide + x]);
    } else { /* big_endian */
      tmp=pixelbuf[y * wide + x];
    } /* swap_byte_order or not */

      pixels[yvalr + x] = tmp >> 22;
      pixels[yvalg + x] = (tmp >> 12) & 1023;
      pixels[yvalb + x] = (tmp >>  2) & 1023;

  } /* x */
 } /* y */

} /* dpx_read_short */

/***********************************************************************************************************/
void
dpx_read (char *inname, float **pixels_read, short *width, short *height, short cineon)
{

FILE *fp;
float *pixels, *pxls_copy;
short x,y;
float red, grn, blu;
short rr, gg, bb;
unsigned int *pixelbuf;
float *pxlbuf;
unsigned int *pp;
unsigned int tmp, hdr_offset;
short wide, tall, ww, hh;
int cinhdr[64];
char dpx_header[2048];
short swap_byte_order;
int yvalr, yvalg, yvalb, yval1, yval2, yval3, yv1, yv2, yv3;
short dpx_is_float = 0;

  swap_byte_order = 0;

  if ((fp = fopen(inname, "rb")) == NULL) {
    if (cineon != 0) {
      printf(" Cannot open cineon input file %s, aborting\n", inname);
    } else { /* dpx */
      printf(" Cannot open dpx input file %s, aborting\n", inname);
	} /* cineon vs dpx */
    exit(1);
  }

  if (cineon != 0) { /* cineon file */
  
    fread(cinhdr, 1, 256 /* 64 ints */, fp);

	tmp = cinhdr[0];
	if (tmp == CIN_MAGIC ) {
#ifdef VERBOSE_PROCESSING
      printf(" reading cineon file big endian header \n ");
#endif /* VERBOSE_PROCESSING */
      swap_byte_order = 0;
    } else {
	  if (tmp == CIN_MAGIC_LE ) {
#ifdef VERBOSE_PROCESSING
        printf(" reading cineon file little endian header \n ");
#endif /* VERBOSE_PROCESSING */
        swap_byte_order = 1;
	  } else { /* magic no good */
	    printf(" cineon magic number = %x no good, aborting\n", tmp);
		exit(1);
	  } /* cin_magic_le or not */
    } /* cin_magic (big endian) or not */


    if (swap_byte_order == 1) {
      hdr_offset = INT_SW(cinhdr[1]);
      tmp  = INT_SW(cinhdr[50]);
      wide = tmp; /* uint to short */
      tmp  = INT_SW(cinhdr[51]);
      tall = tmp; /* uint to short */
    } else { /* big_endian */
      hdr_offset = cinhdr[1];
      wide  = cinhdr[50]; /* uint to short */
      tall  = cinhdr[51]; /* uint to short */
    } /* swap_byte_order or not */

  } else { /* dpx file */

  fread( dpx_header, 1, 2048, fp); 

  if ( *((int *) &dpx_header[0]) == 0x53445058) { /* big-endian "SDPX" in ascii */
#ifdef VERBOSE_PROCESSING
    printf(" reading dpx big-endian file header \n ");
#endif /* VERBOSE_PROCESSING */
    swap_byte_order = 0;
  } else { /* not big-endian header */
    if ( *((int *) &dpx_header[0]) == 0x58504453) { /* little-endian "SDPX" in ascii */
#ifdef VERBOSE_PROCESSING
      printf(" reading dpx little-endian file header \n ");
#endif /* VERBOSE_PROCESSING */
      swap_byte_order = 1;
    } else {
	  printf(" bad magic number = %x in dpx header read, aborting\n", *((int *) &dpx_header[0]));
      exit(1);
    } /* little-endian header or not */
  } /* big-endian header or not */

  if (swap_byte_order == 1) {
    tmp  = INT_SW(*((int *) &dpx_header[772]));
    wide = tmp; /* uint to short */
    tmp  = INT_SW(*((int *) &dpx_header[776]));
    tall = tmp; /* uint to short */
  } else { /* big_endian */
    tmp = *((int *) &dpx_header[772]);
    wide = tmp; /* uint to short */
    tmp = *((int *) &dpx_header[776]);
    tall = tmp; /* uint to short */
  } /* swap_byte_order or not */

  if (dpx_header[803] == 10) {
#ifdef VERBOSE_DECODE
    printf(" dpx packing is normal 10-bit \n ");
#endif /* VERBOSE_DECODE */
  } else {
    if (dpx_header[803] == 32) {
      dpx_is_float = 1;
#ifdef VERBOSE_DECODE
      printf(" dpx packing is 32-bit float \n ");
#endif /* VERBOSE_DECODE */
    } else {
      if (dpx_header[803] == 16) {
        printf(" dpx packing of file %s is 16-bit, which is not (yet) supported (although it should be), aborting\n", inname);
        exit(1);
      } else { /* not 16bits */
        if (dpx_header[803] == 16) {
          printf(" dpx packing of file %s is 12-bit, which is not (yet) supported (although it should be), aborting\n", inname);
          exit(1);
        } else { /* not 12bits */
          printf(" dpx packing of file %s is %d-bits, which is not supported, aborting\n", inname, dpx_header[803]);
          exit(1);
        } /* 12bits or not */
      } /* 16bits or not */
    } /* 32bit float packing or not */
  } /* 10-bit packing or not */


  if (swap_byte_order == 1) {
    hdr_offset = INT_SW(*((unsigned int *) &dpx_header[4]));
  } else { /* big_endian */
    hdr_offset = *((unsigned int *) &dpx_header[4]);
  } /* swap_byte_order or not */

 } /* cineon vs. dpx header */

  printf(" reading file %s width = %d, height = %d, cineon = %d\n", inname, wide, tall, cineon);
  *width = wide;
  *height = tall;

  fseek( fp, hdr_offset, SEEK_SET);

#ifdef VERBOSE_PROCESSING 
  printf(" allocating memory for dpx or cineon pixels \n");
#endif /* VERBOSE_PROCESSING */

 if (dpx_is_float == 1) {
  pxlbuf = (float *) malloc(wide*tall*12); /* rgb float */
 } else { /* not float, must be 10bit */
  pixelbuf = (unsigned int *) malloc(wide*tall*4); /* 4 bytes per pixel in 10-bit dpx packing */
 } /* float or not */

  pixels = (float *) malloc(wide * tall * 12); /* 4-bytes/float * 3-colors */
  *pixels_read = &pixels[0];

#ifdef VERBOSE_PROCESSING 
 printf(" reading dpx or cineon file pixels \n ");
#endif /* VERBOSE_PROCESSING */

  if (dpx_is_float == 1) {
    fread( pxlbuf, 1, wide * tall * sizeof(float) * 3, fp);
    if (ferror(fp)) {
      printf(" error reading data for file %s in dpx_read float, aborting\n", y, inname);
      exit(1);
    }
  } else { /* not float, must be 10bit */
    for( y = 0; y < tall; y++) {
      pp = &pixelbuf[0] + y * wide;
      fread( pp, 1, wide * sizeof(unsigned int), fp);
      if (ferror(fp)) {
        printf(" error reading data for scanline %d for file %s in dpx_read (10bit dpx file), aborting\n", y, inname);
        exit(1);
      }
    } /* y */
  } /* dpx_is_float or not */

 fclose(fp);

  if (dpx_is_float == 1) {

 for( y = 0; y < tall; y++) {
  yvalr = (0*tall + y)*wide;
  yvalg = (1*tall + y)*wide;
  yvalb = (2*tall + y)*wide;
  for( x = 0; x < wide ; x++) {

    if (swap_byte_order == 1) {
unsigned int ttmp;
      tmp=*((unsigned int *) &pxlbuf[(y * wide + x) * 3]);
      ttmp = INT_SW(tmp);
      red = *((float *) &ttmp);
      tmp=*((unsigned int *) &pxlbuf[(y * wide + x) * 3 + 1]);
      ttmp = INT_SW(tmp);
      grn = *((float *) &ttmp);
      tmp=*((unsigned int *) &pxlbuf[(y * wide + x) * 3 + 2]);
      ttmp = INT_SW(tmp);
      blu = *((float *) &ttmp);
    } else { /* dont swap byte order */
      red=pxlbuf[(y * wide + x) * 3];
      grn=pxlbuf[(y * wide + x) * 3 + 1];
      blu=pxlbuf[(y * wide + x) * 3 + 2];
    } /* swap_byte_order or not */

        pixels[yvalr + x] = red;
        pixels[yvalg + x] = grn;
        pixels[yvalb + x] = blu;

  } /* x */
 } /* y */

  } else { /* not float, must be 10bit */

 for( y = 0; y < tall; y++) {
  yvalr = (0*tall + y)*wide;
  yvalg = (1*tall + y)*wide;
  yvalb = (2*tall + y)*wide;
  for( x = 0; x < wide ; x++) {

    if (swap_byte_order == 1) {
      tmp=INT_SW(pixelbuf[y * wide + x]);
    } else { /* big_endian */
      tmp=pixelbuf[y * wide + x];
    } /* swap_byte_order or not */

      rr = tmp >> 22;
      gg = (tmp >> 12) & 1023;
      bb = (tmp >>  2) & 1023;

      red = rr/1023.0;
      grn = gg/1023.0;
      blu = bb/1023.0;

      pixels[yvalr + x] = red;
      pixels[yvalg + x] = grn;
      pixels[yvalb + x] = blu;

  } /* x */
 } /* y */
 } /* dpx_is_float or not */

  if (dpx_is_float == 1) {
    free(pxlbuf);
   } else { /* not float, must be 10bit */
    free(pixelbuf);
  } /* dpx_is_float or not */

} /* dpx_read */

#undef INTEL_LE /* it is arbitrary which endian the dpx file is written, since the magic in the header tells the reader which to read (and must support both) */
#define HEADER_SIZE 2048 /* can also make this 8192, or other */

/***********************************************************************************************************/
void
dpx_write_10bit_from_float(char *outname, float *pixel_result, short width, short height)
{

FILE *fp_out;

unsigned int tmp, size;
short ired, igrn, iblu, x, y, i;
float red,grn,blu;

unsigned int *pixelbuf;
unsigned char dpx_header[HEADER_SIZE];
int yvalr, yvalg, yvalb, yval0;


     if ((fp_out = fopen(outname,"wb" /* rb for reading, wb for writing */)) == NULL) {
        printf(" Cannot open dpx output file %s \n", outname);
        exit(1);
      } /* fopen */

  for (i=0; i< HEADER_SIZE; i++) { /* clear header array */
    dpx_header[i]=0;
  } /* i */

  pixelbuf = (unsigned int *) malloc(width*height*4); /* 4 bytes per pixel in 10-bit dpx packing */


 for( y = 0; y < height; y++) {
  yval0 = y * width;
  yvalr = (0*height + y)*width;
  yvalg = (1*height + y)*width;
  yvalb = (2*height + y)*width;
  for( x = 0; x < width ; x++) {

      red = MAX(0.0, MIN(1023.0, 1023.0 * pixel_result[yvalr + x]));
      grn = MAX(0.0, MIN(1023.0, 1023.0 * pixel_result[yvalg + x]));
      blu = MAX(0.0, MIN(1023.0, 1023.0 * pixel_result[yvalb + x]));

	  ired = red; /* float to short */
      igrn = grn;
      iblu = blu;

      tmp = (ired << 22) | (igrn << 12) | (iblu << 2);

#ifdef INTEL_LE
      pixelbuf[yval0 + x] = INT_SW(tmp);
#else /* not INTEL_LE */
      pixelbuf[yval0 + x] = tmp;
#endif /* INTEL_LE or not */

  } /* x */
 } /* y */



#ifdef VERBOSE_PROCESSING
printf(" writing dpx file  %s having width = %d height = %d\n", outname, width, height);
#endif /* VERBOSE_PROCESSING */


#ifdef INTEL_LE
  tmp = width; /* short to uint */
  *((int *) &dpx_header[772])   = INT_SW(tmp);
  tmp = height; /* short to uint */
  *((int *) &dpx_header[776]) = INT_SW(tmp);
  tmp = width; /* short to uint */
  *((int *) &dpx_header[1424]) = INT_SW(tmp);
  tmp = height; /* short to uint */
  *((int *) &dpx_header[1428]) = INT_SW(tmp);
#else /* not INTEL_LE */
  *((int *) &dpx_header[772]) = width;
  *((int *) &dpx_header[776]) = height;
  *((int *) &dpx_header[1424]) = width;
  *((int *) &dpx_header[1428]) = height;
#endif /* INTEL_LE or not */

  size = width * height * 4; /* 4 bytes per pixel, (in bytes???) is this right for the size field???, or should this be in pixels??? */

#ifdef INTEL_LE
  tmp = 2048/*8192*/; /* 6144 is the pad since the header is 2048 (pad needed for photoshop, but not needed for graphic_converter) */
  *((int *) &dpx_header[4]) = INT_SW(tmp); /* image offset */
  tmp = size + tmp;
  *((int *) &dpx_header[16]) = INT_SW(tmp); /* total file size */
  tmp = 0x53445058;  /* = 1396985944 decimal, write big-endian file, even if little-endian/Intel */
  *((int *) &dpx_header[0]) = INT_SW(tmp);
  sprintf(((char *) &dpx_header[8]), "V1.0    "); /* version, not sure about little-endian byte order for this character string */

  *((int *) &dpx_header[20]) = 0; /* ditto key */
  tmp = 1664;
  *((int *) &dpx_header[24]) = INT_SW(tmp); /* generic length */
  tmp = 384;
  *((int *) &dpx_header[28]) = INT_SW(tmp); /* industry specific length */
  *((int *) &dpx_header[32]) = 0; /* user length */
  *((int *) &dpx_header[660]) = 0xffffffff; /* encryption key */
#else /* not INTEL_LE */
  *((int *) &dpx_header[4]) = 2048/*8192*/; /* image offset, 6144 is the pad since the header is 2048 (pad needed for photoshop, but not needed for graphic_converter) */
  *((int *) &dpx_header[16]) = size + 2048; /* total file size */
  *((int *) &dpx_header[0]) = 0x53445058  /* = 1396985944 decimal */;
  sprintf(((char *) &dpx_header[8]), "V1.0    "); /* version */
  *((int *) &dpx_header[20]) = 0; /* ditto key */
  *((int *) &dpx_header[24]) = 1664; /* generic length */
  *((int *) &dpx_header[28]) = 384; /* industry specific length */
  *((int *) &dpx_header[32]) = 0; /* user length */
  *((int *) &dpx_header[660]) = 0xffffffff; /* encryption key */
#endif /* INTEL_LE or not */



#ifdef INTEL_LE
 *((short *) &dpx_header[768]) = 0; /* orientation */
 tmp = 1;
 *((short *) &dpx_header[770]) = SHORT_SW(tmp); /* number of elements */

 *((int  *) &dpx_header[780]) = 0; /* data signed */
 dpx_header[800] = 50; /* descriptor, rgb data */
 dpx_header[801] = 6; /* transfer characteristic, 1 is printing density, 6 is Rec709 */
 dpx_header[802] = 6; /* colorimetric, 1 is printing density, 6 is Rec709 */
 tmp = 10;
 dpx_header[803] = 10; /* bits per element */
 *((short *) &dpx_header[804]) = 1; /* packing, packed into 32-bit words */
 *((short *) &dpx_header[806]) = 0; /* encoding, no run-length */
 tmp = 2048/*8192*/; /* byte offset to red pixels */
 *((int *) &dpx_header[808]) = INT_SW(tmp); /* offset to image */
 *((int *) &dpx_header[812]) = 0; /* no end of line padding */
 *((int *) &dpx_header[816]) = 0; /* no end of image padding */
#else /* not INTEL_LE */
 *((short *) &dpx_header[768]) = 0; /* orientation */
 *((short *) &dpx_header[770]) = 1; /* number of elements */

 *((int  *) &dpx_header[780]) = 0; /* data signed */
 dpx_header[800] = 50; /* descriptor, rgb data */
 dpx_header[801] = 6; /* transfer characteristic, 1 is printing density, 6 is Rec709 */
 dpx_header[802] = 6; /* colorimetric, 1 is printing density, 6 is Rec709 */
 dpx_header[803] = 10; /* bits per element */
 *((short *) &dpx_header[804]) = 0; /* packing, packed into 32-bit words */
 *((short *) &dpx_header[806]) = 0; /* encoding, no run-length */
 *((int *) &dpx_header[808]) = 2048/*8192*/; /* byte offset to red pixels */
 *((int *) &dpx_header[812]) = 0; /* no end of line padding */
 *((int *) &dpx_header[816]) = 0; /* no end of image padding */
#endif /* INTEL_LE or not */

#ifdef VERBOSE_PROCESSING
    printf(" writing dpx file header of size 2048 \n ");
#endif /* VERBOSE_PROCESSING */

    fwrite(dpx_header, 1, 2048, fp_out); 

//    fwrite(pixelbuf/* dummy */,1,8192 - 2048,fp_out); /* pad (needed for photoshop, but not needed for graphic_converter) */

#ifdef VERBOSE_PROCESSING
    printf(" writing dpx pixels \n ");
#endif /* VERBOSE_PROCESSING */

    fwrite( pixelbuf, 1, height*width*4, fp_out);

#ifdef VERBOSE_PROCESSING
    printf(" closing dpx output file \n");
#endif /* VERBOSE_PROCESSING */

    fclose(fp_out);

    free(pixelbuf);

} /* write_dpx */
/*****************************************************************************************************************/

void
dpx_write_float(char *outname, float *pixel_result, short width, short height)
{

FILE *fp_out;

unsigned int tmp, size;
short x, y, i;
float red, grn, blu;

float *pixelbuf;
unsigned char dpx_header[HEADER_SIZE];
int yvalr, yvalg, yvalb, yval0;


     if ((fp_out = fopen(outname,"wb" /* rb for reading, wb for writing */)) == NULL) {
        printf(" Cannot open dpx output file %s \n", outname);
        exit(1);
      } /* fopen */

  for (i=0; i< HEADER_SIZE; i++) { /* clear header array */
    dpx_header[i]=0;
  } /* i */


  pixelbuf = (float *) malloc(width*height*12); /* rgb floats */


 for( y = 0; y < height; y++) {
  yval0 = y * width * 3;
  yvalr = (0*height + y)*width;
  yvalg = (1*height + y)*width;
  yvalb = (2*height + y)*width;
  for( x = 0; x < width ; x++) {

      red = pixel_result[yvalr + x];
      grn = pixel_result[yvalg + x];
      blu = pixel_result[yvalb + x];

#ifdef INTEL_LE
unsigned int ttmp, ttmp2;
      ttmp = *((unsigned int *) &red);
      ttmp2 = INT_SW(ttmp);
      red = *((float *) &ttmp2);

      ttmp = *((unsigned int *) &grn);
      ttmp2 = INT_SW(ttmp);
      grn = *((float *) &ttmp2);

      ttmp = *((unsigned int *) &blu);
      ttmp2 = INT_SW(ttmp);
      blu = *((float *) &ttmp2);

#endif /* INTEL_LE or not */

      pixelbuf[yval0 + 3 * x    ] = red;
      pixelbuf[yval0 + 3 * x + 1] = grn;
      pixelbuf[yval0 + 3 * x + 2] = blu;

  } /* x */
 } /* y */



#ifdef VERBOSE_DECODE
printf(" writing dpx file  %s having width = %d height = %d\n", outname, width, height);
#endif /* VERBOSE_DECODE */

/* fill undefined values */
  for(i=784; i<800;   i++) { dpx_header[i]=0xff; }
  for(i=812; i<816;   i++) { dpx_header[i]=0xff; }
  for(i=852; i<892;   i++) { dpx_header[i]=0xff; }
  for(i=924; i<964;   i++) { dpx_header[i]=0xff; }
  for(i=996; i<1036;  i++) { dpx_header[i]=0xff; }
  for(i=1068; i<1108; i++) { dpx_header[i]=0xff; }
  for(i=1140; i<1180; i++) { dpx_header[i]=0xff; }
  for(i=1212; i<1252; i++) { dpx_header[i]=0xff; }
  for(i=1284; i<1324; i++) { dpx_header[i]=0xff; }
  for(i=1408; i<1432; i++) { dpx_header[i]=0xff; }
  for(i=1620; i<1644; i++) { dpx_header[i]=0xff; }
  for(i=1712; i<1732; i++) { dpx_header[i]=0xff; }
  for(i=1920; i<1972; i++) { dpx_header[i]=0xff; }
  dpx_header[1931]=0; /* this one is not undefined, it defines "0" (zero) for byte alignment */

#ifdef INTEL_LE
  tmp = width; /* short to uint */
  *((int *) &dpx_header[772])   = INT_SW(tmp);
  tmp = height; /* short to uint */
  *((int *) &dpx_header[776]) = INT_SW(tmp);
  tmp = 0xffffffff /* width */; /* short to uint */
  *((int *) &dpx_header[1424]) = INT_SW(tmp);
  tmp = 0xffffffff /* height */; /* short to uint */
  *((int *) &dpx_header[1428]) = INT_SW(tmp);
#else /* not INTEL_LE */
  *((int *) &dpx_header[772]) = width;
  *((int *) &dpx_header[776]) = height;
  *((int *) &dpx_header[1424]) = 0xffffffff /* width */;
  *((int *) &dpx_header[1428]) = 0xffffffff /* height */;
#endif /* INTEL_LE or not */

  size = width * height * 12; /* 12 bytes per pixel, (in bytes???) is this right for the size field???, or should this be in pixels??? */

#ifdef INTEL_LE
  tmp = 0x53445058;  /* = 1396985944 decimal, write big-endian file, even if little-endian/Intel, although turning off intel_le on intel machines will write little_endian */
  *((int *) &dpx_header[0]) = INT_SW(tmp);
  sprintf(((char *) &dpx_header[8]), "v2.0"); /* version, not sure about little-endian byte order for this character string */
  tmp = HEADER_SIZE;
  *((int *) &dpx_header[4]) = INT_SW(tmp); /* image offset */
  tmp = size + tmp;
  *((int *) &dpx_header[16]) = INT_SW(tmp); /* total file size */

//  *((int *) &dpx_header[20]) = 0; /* ditto key */
//  tmp = 1664;
//  *((int *) &dpx_header[24]) = INT_SW(tmp); /* generic length */
//  tmp = 384;
//  *((int *) &dpx_header[28]) = INT_SW(tmp); /* industry specific length */
//  *((int *) &dpx_header[32]) = 0; /* user length */
  for(i=20; i<36; i++) { dpx_header[i]=0xff; /* make ditto key, generic length, industry specific length, and user length all be undefined */ }

  *((int *) &dpx_header[660]) = 0xffffffff; /* encryption key */
#else /* not INTEL_LE */
  *((int *) &dpx_header[0]) = 0x53445058  /* = 1396985944 decimal */;
  sprintf(((char *) &dpx_header[8]), "v2.0"); /* version */
  *((int *) &dpx_header[4]) = HEADER_SIZE;
  *((int *) &dpx_header[16]) = size + HEADER_SIZE; /* total file size */
//  *((int *) &dpx_header[20]) = 0; /* ditto key */
//  *((int *) &dpx_header[24]) = 1664; /* generic length */
//  *((int *) &dpx_header[28]) = 384; /* industry specific length */
//  *((int *) &dpx_header[32]) = 0; /* user length */
  for(i=20; i<36; i++) { dpx_header[i]=0xff; /* make ditto key, generic length, industry specific length, and user length all be undefined */ }

  *((int *) &dpx_header[660]) = 0xffffffff; /* encryption key */
#endif /* INTEL_LE or not */



#ifdef INTEL_LE
 *((short *) &dpx_header[768]) = 0; /* orientation */
 tmp = 1;
 *((short *) &dpx_header[770]) = SHORT_SW(tmp); /* number of elements */

 *((int  *) &dpx_header[780]) = 1; /* data signed */
 dpx_header[800] = 50; /* descriptor, rgb data */
 dpx_header[801] = 6; /* transfer characteristic, 6 = video_gamma, 1 = printing density */
 dpx_header[802] = 6; /* colorimetric, 6 = video_gamma, 1 = printing density */
 tmp = 10;
 dpx_header[803] = 32; /* bits per element */
 *((short *) &dpx_header[804]) = 0; /* packing, no packing (into 32-bit floats) */
 *((short *) &dpx_header[806]) = 0; /* encoding, no run-length */
 tmp = HEADER_SIZE; /* byte offset to red pixels */
 *((int *) &dpx_header[808]) = INT_SW(tmp); /* offset to image */
 *((int *) &dpx_header[812]) = 0; /* no end of line padding */
 *((int *) &dpx_header[816]) = 0; /* no end of image padding */
#else /* not INTEL_LE */
 *((short *) &dpx_header[768]) = 0; /* orientation */
 *((short *) &dpx_header[770]) = 1; /* number of elements */

 *((int  *) &dpx_header[780]) = 1; /* data signed */
 dpx_header[800] = 50; /* descriptor, rgb data */
 dpx_header[801] = 6; /* transfer characteristic, 6 = video_gamma, 1 = printing density */
 dpx_header[802] = 6; /* colorimetric, 6 = video_gamma, 1 = printing density */
 dpx_header[803] = 32; /* bits per element */
 *((short *) &dpx_header[804]) = 0; /* packing, no packing (into 32-bit floats) */
 *((short *) &dpx_header[806]) = 0; /* encoding, no run-length */
 *((int *) &dpx_header[808]) = HEADER_SIZE; /* byte offset to red pixels */
 *((int *) &dpx_header[812]) = 0; /* no end of line padding */
 *((int *) &dpx_header[816]) = 0; /* no end of image padding */
#endif /* INTEL_LE or not */

#ifdef VERBOSE_DECODE
    printf(" writing dpx file header of size %d \n ", HEADER_SIZE);
#endif /* VERBOSE_DECODE */

    fwrite(dpx_header, 1, HEADER_SIZE, fp_out); 

    if (ferror(fp_out)) {
      printf(" error writing header for file %s in dpx_write, aborting\n", outname);
      exit(1);
    }

#ifdef VERBOSE_DECODE
    printf(" writing dpx pixels \n ");
#endif /* VERBOSE_DECODE */

    fwrite( pixelbuf, 1, height*width*12, fp_out);

    if (ferror(fp_out)) {
      printf(" error writing data to file %s in dpx_write, aborting\n", outname);
      exit(1);
    }

#ifdef VERBOSE_DECODE
    printf(" closing dpx32 float output file \n");
#endif /* VERBOSE_DECODE */

 if(fclose(fp_out)) {
   printf(" error closing file %s in dpx_write_float, aborting\n", outname);
   exit(1);
 }

    free(pixelbuf);


} /* dpx_write_float */


/***********************************************************************************************************/
int
main(int argc, char **argv)
{
short odt_type=0;
short x,y,c,i;
char outfile[300], infile[300];
short num_chars;
int first, last, frame;

short h_reso, v_reso, ans_ready;
short open_exr_vs_dpx_out = 0; /* 0 is open_exr, 1 is dpx */
float s, t, tmp;
int ii;
float rOut, gOut, bOut;
float rgbOut[3], aces[3];


if (argc < 6) {
  printf(" usage: %s infiles(dpx), outfiles(exr), first_frame, last_frame, device_type(1=GD10_Rec709_MDR, 2=GD10_p3_d60_HDR)\n", argv[0]);
  exit(1);
}

 first = atoi(argv[3]);
 last  = atoi(argv[4]);

 odt_type = atoi(argv[5]);

 printf(" inverse (dpx to exr) processing frames %d to %d with device_type = %d\n", first, last, odt_type);

/* warning, 10bits for the lookup dimension gets an overflow in these mallocs, additional work would be needed for 10bit mallocs */
    count     = (short *) malloc(LOOKUP_DIMENSION * LOOKUP_DIMENSION * LOOKUP_DIMENSION * 2); /* 2bytes/short */
    lookup[0] = (float *) malloc(LOOKUP_DIMENSION * LOOKUP_DIMENSION * LOOKUP_DIMENSION * 4); /* 4bytes/float */
    lookup[1] = (float *) malloc(LOOKUP_DIMENSION * LOOKUP_DIMENSION * LOOKUP_DIMENSION * 4); /* 4bytes/float */
    lookup[2] = (float *) malloc(LOOKUP_DIMENSION * LOOKUP_DIMENSION * LOOKUP_DIMENSION * 4); /* 4bytes/float */


/* note: multi-threading the following is not easy, since atomic operations are required (and usually not available in multi-threading tools) */
{
 unsigned int touch_count = 0, neighbor_touch_count = 0;
 short i,j,k, ii,jj,kk;
 for(i=0; i<LOOKUP_DIMENSION; i++) {
  for(j=0; j<LOOKUP_DIMENSION; j++) {
   for(k=0; k<LOOKUP_DIMENSION; k++) {
     count[i*LOOKUP_DIMENSION*LOOKUP_DIMENSION + j*LOOKUP_DIMENSION + k]=0; /* initialize as empty */
   } /* k */
  } /* j */
 } /* i */


 float random_batch[RANDOM_BATCH_SIZE][3];
 for(jj=0; jj<NUMBER_OF_RANDOM_BATCHES; jj++) {
  float range_top_r, range_top_g, range_top_b;
  if ((jj & 7) == 0) { /* alternate between different ranges */
   range_top_r = range_top_g = range_top_b = RANGE_TOP_A;
  } else {
   if ((jj & 7) == 1) { /* alternate between different ranges */
    range_top_r = range_top_g = range_top_b = RANGE_TOP_B;
   } else {
    if ((jj & 7) == 2) { /* alternate between different ranges */
     range_top_r = range_top_g = range_top_b = RANGE_TOP_C;
    } else {
     if ((jj & 7) == 3) { /* alternate between different ranges */
      range_top_r = range_top_g = range_top_b = RANGE_TOP_D;
     } else {
      if ((jj & 7) == 4) {
       range_top_r = RNG_RR;
       range_top_g = RNG_GR;
       range_top_b = RNG_BR;
      } else {
       if ((jj & 7) == 5) {
        range_top_r = RNG_RG;
        range_top_g = RNG_GG;
        range_top_b = RNG_BG;
       } else {
        if ((jj & 7) == 6) {
         range_top_r = RNG_RB;
         range_top_g = RNG_GB;
         range_top_b = RNG_BB;
        } else { /* 7 */
         range_top_r = RNG_LG;
         range_top_g = RNG_LG;
         range_top_b = RNG_LG;
        } /* 6 or 7 */
       } /* 5 or not */
      } /* 4 or not */
     } /* 3 or not */
    } /* 2 or not */
   } /* 1 or not */
  } /* 0 or not */
  for(i=0; i<RANDOM_BATCH_SIZE; i++) {
   random_batch[i][0]= (range_top_r*rand())/(1.0*RAND_MAX);
   random_batch[i][1]= (range_top_g*rand())/(1.0*RAND_MAX);
   random_batch[i][2]= (range_top_b*rand())/(1.0*RAND_MAX);
  } /* i */

  for(i=0; i<RANDOM_BATCH_SIZE; i++) {
   for(j=0; j<RANDOM_BATCH_SIZE; j++) {
    for(k=0; k<RANDOM_BATCH_SIZE; k++) {

float rr, gg, bb;
      aces[0] = random_batch[i][0];
      aces[1] = random_batch[j][1];
      aces[2] = random_batch[k][2];

      rr = aces[0]; gg = aces[1]; bb = aces[2];

      RRT_GD10(aces, rgbOut)

      if (odt_type == 1) {

        process_odt_gd9_Rec709_g2pt4_MDR

      } else { /* odt_type != 1 */
        if (odt_type == 2) {

          process_odt_gd9_p3_d60_g2pt4_HDR

        } else { /* odt_type != 2 */
          printf(" unknown odt_type = %d (1=GD10_Rec709_MDR, 2=GD10_p3_d60_HDR), aborting\n", odt_type);
          exit(1);
        } /* 2 or not */
      } /* 1 or not */

short r_int, g_int, b_int;
      r_int = rOut*1023.499; /* float to short int */
      g_int = gOut*1023.499;
      b_int = bOut*1023.499;
if((r_int > 1023) || (g_int > 1023) || (b_int > 1023)) {printf(" error, out of 1023 range r_int=%d g_int=%d b_int=%d\n",r_int, g_int, b_int);exit(1);}

#if 1
{
        r_int = (MIN(1023, (r_int+2)) >> 2); /* 8bit from 10bit with crude rounding */
        g_int = (MIN(1023, (g_int+2)) >> 2); /* 8bit from 10bit with crude rounding */
        b_int = (MIN(1023, (b_int+2)) >> 2); /* 8bit from 10bit with crude rounding */
#else
      if (((r_int & 3)==0) && ((g_int & 3)==0) && ((b_int & 3)==0)) { /* ignore 63 of 64 pixels that dont fall exactly on 10bit points for 8bit lookup/interpolate */
        r_int = (r_int >> 2); /* 8bit from 10bit */
        g_int = (g_int >> 2); /* 8bit from 10bit */
        b_int = (b_int >> 2); /* 8bit from 10bit */
#endif
int index;
        index = r_int * LOOKUP_DIMENSION * LOOKUP_DIMENSION + g_int * LOOKUP_DIMENSION + b_int;
        if (count[index] <= 0) { /* first direct touch of this index */
         if (count[index] < 0) { neighbor_touch_count--; /*remove from neighbor touch count, and next replace with direct touch count */ }
         touch_count++; /* count all indices touched to indicate relative progress */
         count[index]=1;
         lookup[0][index] = rr;
         lookup[1][index] = gg;
         lookup[2][index] = bb;
        } else { /* not first */
         float lower_count = count[index];
         count[index]++;
         float higher_count = count[index];
/* average in new value */
         lookup[0][index] = (lookup[0][index] * lower_count) / higher_count + rr / higher_count;
         lookup[1][index] = (lookup[1][index] * lower_count) / higher_count + gg / higher_count;
         lookup[2][index] = (lookup[2][index] * lower_count) / higher_count + bb / higher_count;
        }
/* touch all nearby values not having any count, in case they are never touched, but dont update count */
short r_int_minus1, g_int_minus1, b_int_minus1, r_int_plus1, g_int_plus1, b_int_plus1;
        r_int_minus1 = MAX(0, r_int-1);
        g_int_minus1 = MAX(0, g_int-1);
        b_int_minus1 = MAX(0, b_int-1);
        r_int_plus1 = MIN(255, r_int+1);
        g_int_plus1 = MIN(255, g_int+1);
        b_int_plus1 = MIN(255, b_int+1);

#define touch_neighbors(rint, gint, bint) \
        index = rint * LOOKUP_DIMENSION * LOOKUP_DIMENSION + gint * LOOKUP_DIMENSION + bint; \
        if (count[index] <= 0) { /* no direct touches of this index yet */ \
          if (count[index]==0) { /* first neigbhor touch of this index */ \
            neighbor_touch_count++; /* count all neighbor indices touched */ \
            count[index]=-1; /* keep negative count of all of these neightbor touches */ \
            lookup[0][index] = rr; \
            lookup[1][index] = gg; \
            lookup[2][index] = bb; \
          } else { /* not first */ \
            float lower_count = -count[index]; \
            count[index]--; \
            float higher_count = -count[index]; \
/* average in new value */ \
            lookup[0][index] = (lookup[0][index] * lower_count) / higher_count + rr / higher_count; \
            lookup[1][index] = (lookup[1][index] * lower_count) / higher_count + gg / higher_count; \
            lookup[2][index] = (lookup[2][index] * lower_count) / higher_count + bb / higher_count; \
          } /* count[index]==0 or not */\
        } /* count[index] <= 0 */

        touch_neighbors(r_int, g_int,        b_int_minus1)
        touch_neighbors(r_int, g_int,        b_int_plus1)
        touch_neighbors(r_int, g_int_minus1, b_int)
        touch_neighbors(r_int, g_int_plus1,  b_int)
        touch_neighbors(r_int, g_int_minus1, b_int_minus1)
        touch_neighbors(r_int, g_int_plus1,  b_int_minus1)
        touch_neighbors(r_int, g_int_minus1, b_int_plus1)
        touch_neighbors(r_int, g_int_plus1,  b_int_plus1)

        touch_neighbors(r_int_minus1, g_int,        b_int)
        touch_neighbors(r_int_minus1, g_int,        b_int_minus1)
        touch_neighbors(r_int_minus1, g_int,        b_int_plus1)
        touch_neighbors(r_int_minus1, g_int_minus1, b_int)
        touch_neighbors(r_int_minus1, g_int_plus1,  b_int)
        touch_neighbors(r_int_minus1, g_int_minus1, b_int_minus1)
        touch_neighbors(r_int_minus1, g_int_plus1,  b_int_minus1)
        touch_neighbors(r_int_minus1, g_int_minus1, b_int_plus1)
        touch_neighbors(r_int_minus1, g_int_plus1,  b_int_plus1)

        touch_neighbors(r_int_plus1,  g_int, b_int)
        touch_neighbors(r_int_plus1,  g_int, b_int_minus1)
        touch_neighbors(r_int_plus1,  g_int, b_int_plus1)
        touch_neighbors(r_int_plus1,  g_int_minus1, b_int)
        touch_neighbors(r_int_plus1,  g_int_plus1,  b_int)
        touch_neighbors(r_int_plus1,  g_int_minus1, b_int_minus1)
        touch_neighbors(r_int_plus1,  g_int_plus1,  b_int_minus1)
        touch_neighbors(r_int_plus1,  g_int_minus1, b_int_plus1)
        touch_neighbors(r_int_plus1,  g_int_plus1,  b_int_plus1)


      } /* low 2bits all zero */      
   } /* k */
  } /* j */
 } /* i */
   printf(" finished batch %d of %d, have now touched %d, and touched neighbors %d (total %d) out of %d\n",
     jj+1, NUMBER_OF_RANDOM_BATCHES, touch_count, neighbor_touch_count, touch_count + neighbor_touch_count, LOOKUP_DIMENSION * LOOKUP_DIMENSION * LOOKUP_DIMENSION);
 } /* jj */

#define TOP_START_VAL 5.0

/* fill in the holes */
 {
   short r,g,b;
   float first_red, first_grn, first_blu, red, grn, blu;
int val_indx;

/* note that this order of loops is not optimal */
/* fill the holes */
int fill_count = 0, valid_count = 0;
#define HOLE_ITERATIONS 7 /* number of iterations starting from the bottom (black) and sweeping up */
  red = first_red = 0.0;
  grn = first_grn = 0.0;
  blu = first_blu = 0.0;
  for(i=0; i<HOLE_ITERATIONS; i++) {
   short top_amt = (LOOKUP_DIMENSION>>(HOLE_ITERATIONS-i-1));
   for (r=0; r<top_amt; r++) {
    for (g=0; g<top_amt; g++) {
     for (b=0; b<top_amt; b++) {
      val_indx = r * LOOKUP_DIMENSION * LOOKUP_DIMENSION + g * LOOKUP_DIMENSION + b;
      if (count[val_indx] != 0) { /* this doesnt distinguish between positive, which is direct touch, vs negative neighbor touch */
       valid_count++;
       red = lookup[0][val_indx]; /* update for next hole that appears */
       grn = lookup[1][val_indx]; /* update for next hole that appears */
       blu = lookup[2][val_indx]; /* update for next hole that appears */
      } else { /* hole, fill with last */
       fill_count++;
       count[val_indx] = -1; /* mark as filled */
       lookup[0][val_indx] = red; /* fills from the bottom up */
       lookup[1][val_indx] = grn;
       lookup[2][val_indx] = blu;
/* average in any valid neighbors */
short r_int_minus1, g_int_minus1, b_int_minus1, r_int_plus1, g_int_plus1, b_int_plus1;
int index;
        r_int_minus1 = MAX(0, r-1);
        g_int_minus1 = MAX(0, g-1);
        b_int_minus1 = MAX(0, b-1);
        r_int_plus1 = MIN(255, r+1);
        g_int_plus1 = MIN(255, g+1);
        b_int_plus1 = MIN(255, b+1);

#define average_in_valid_neighbors(rr, gg, bb) \
        index = rr * LOOKUP_DIMENSION * LOOKUP_DIMENSION + gg * LOOKUP_DIMENSION + bb; \
        if (count[index] != 0) { /* valid neighbor */ \
          float lower_count = -count[val_indx]; \
          count[val_indx]--; \
          float higher_count = -count[val_indx]; \
/* average in new value */ \
          lookup[0][val_indx] = (lookup[0][val_indx] * lower_count) / higher_count + lookup[0][index] / higher_count; \
          lookup[1][val_indx] = (lookup[1][val_indx] * lower_count) / higher_count + lookup[1][index] / higher_count; \
          lookup[2][val_indx] = (lookup[2][val_indx] * lower_count) / higher_count + lookup[2][index] / higher_count; \
        } /* count[index] > 0 */

        average_in_valid_neighbors(r, g,        b_int_minus1)
        average_in_valid_neighbors(r, g,        b_int_plus1)
        average_in_valid_neighbors(r, g_int_minus1, b)
        average_in_valid_neighbors(r, g_int_plus1,  b)
        average_in_valid_neighbors(r, g_int_minus1, b_int_minus1)
        average_in_valid_neighbors(r, g_int_plus1,  b_int_minus1)
        average_in_valid_neighbors(r, g_int_minus1, b_int_plus1)
        average_in_valid_neighbors(r, g_int_plus1,  b_int_plus1)

        average_in_valid_neighbors(r_int_minus1, g,        b)
        average_in_valid_neighbors(r_int_minus1, g,        b_int_minus1)
        average_in_valid_neighbors(r_int_minus1, g,        b_int_plus1)
        average_in_valid_neighbors(r_int_minus1, g_int_minus1, b)
        average_in_valid_neighbors(r_int_minus1, g_int_plus1,  b)
        average_in_valid_neighbors(r_int_minus1, g_int_minus1, b_int_minus1)
        average_in_valid_neighbors(r_int_minus1, g_int_plus1,  b_int_minus1)
        average_in_valid_neighbors(r_int_minus1, g_int_minus1, b_int_plus1)
        average_in_valid_neighbors(r_int_minus1, g_int_plus1,  b_int_plus1)

        average_in_valid_neighbors(r_int_plus1,  g, b)
        average_in_valid_neighbors(r_int_plus1,  g, b_int_minus1)
        average_in_valid_neighbors(r_int_plus1,  g, b_int_plus1)
        average_in_valid_neighbors(r_int_plus1,  g_int_minus1, b)
        average_in_valid_neighbors(r_int_plus1,  g_int_plus1,  b)
        average_in_valid_neighbors(r_int_plus1,  g_int_minus1, b_int_minus1)
        average_in_valid_neighbors(r_int_plus1,  g_int_plus1,  b_int_minus1)
        average_in_valid_neighbors(r_int_plus1,  g_int_minus1, b_int_plus1)
        average_in_valid_neighbors(r_int_plus1,  g_int_plus1,  b_int_plus1)

        if (count[val_indx] < -1) { /* update running value if valid neighbors */
         red = lookup[0][val_indx]; /* update for next hole that appears */
         grn = lookup[1][val_indx]; /* update for next hole that appears */
         blu = lookup[2][val_indx]; /* update for next hole that appears */
        } /* count[val_indx] < -1 */

      } /* count[val_indx] != 0 vs hole */
     } /* b */
     blu = first_blu;
    } /* g */
    grn = first_grn;
   } /* r */
   red = first_red;
  } /* i */
  printf(" done filling holes, valid count = %d fill_count = %d\n", valid_count, fill_count);
 }
}


/****************************************************************************************************************************************/


   for (frame=first; frame <= last; frame++) {


     sprintf(infile, argv[1], frame);
     num_chars = strlen(infile); /* length of outfile string */
     if ((!strcmp(&infile[num_chars-1], "x"))||(!strcmp(&infile[num_chars-1], "X"))) { /* DPX file ending in ".dpx" */
      printf(" processing input file %s\n", infile);
     } else { /* not dpx */
      printf(" unknown filetype for reading, since extension doesn't end in x, only dpx reading supported, infile = %s, aborting\n", infile);
      exit(1);
     } /* exr or not */

     sprintf(outfile, argv[2], frame);
     num_chars = strlen(outfile); /* length of outfile string */
     if ((!strcmp(&outfile[num_chars-1], "r"))||(!strcmp(&outfile[num_chars-1], "R"))) { /* EXR file ending in ".exr" */
      printf(" processing output file %s\n", outfile);
      open_exr_vs_dpx_out = 0; /* exr */
     } else { /* not exr */
      printf(" unknown filetype for writing, since extension doesn't end in r, only exr writing supported, outfile = %s, aborting\n", outfile);
      exit(1);
     } /* exr or not */

unsigned short *pixels_read = NULL;

printf(" begin reading %s\n", infile);

     dpx_read_short(infile, &pixels_read, &h_reso, &v_reso, 0/*not cineon*/);

 printf(" begin writing %s\n", outfile);

{
    Array2D<Rgba> pixels_out_exr(v_reso, h_reso);

int val_indx;
short r, g, b, rr, gg, bb, r_plus1, g_plus1, b_plus1;
float red, grn, blu, col000[3], col001[3], col010[3], col011[3], col100[3], col101[3], col110[3], col111[3];
float col00[3], col01[3], col10[3], col11[3], col0[3], col1[3];
float red_frac, grn_frac, blu_frac, red_frac_flt, grn_frac_flt, blu_frac_flt;

/* note: the following could benefit from multi-threading */
      for (y=0; y<v_reso; y++) {
        for (x=0; x<h_reso; x++) {
          rr = pixels_read[(0 * v_reso + y) * h_reso + x];
          gg = pixels_read[(1 * v_reso + y) * h_reso + x];
          bb = pixels_read[(2 * v_reso + y) * h_reso + x];

/* limit maximum with hard clip */
#define MAX_LIMIT 1014

/* stop increasing desat amount at these spots (slope change could be handled with sine or other function) */
#define MAX_DESAT_R MIN(MAX_LIMIT, 700.0)
#define MAX_DESAT_G MIN(MAX_LIMIT, 800.0)
#define MAX_DESAT_B MIN(MAX_LIMIT, 800.0)
#define MAX_DESAT_C MIN(MAX_LIMIT, 900.0)
#define MAX_DESAT_M MIN(MAX_LIMIT, 1023.0)
#define MAX_DESAT_Y MIN(MAX_LIMIT, 1023.0)

/* reduce saturation of full-saturation colors */
#define RED_DESAT_AMT 85.0
#define GRN_DESAT_AMT 22.0
#define BLU_DESAT_AMT 5.0
#define CYAN_DESAT_AMT 8.0
#define MAGENTA_DESAT_AMT 0.0
#define YELLOW_DESAT_AMT 0.0

rr = MIN(MAX_LIMIT, rr);
gg = MIN(MAX_LIMIT, gg);
bb = MIN(MAX_LIMIT, bb);

short maxx, minn;
float range, midrange, desat_inc;
         maxx = MAX(rr, MAX(gg, bb));
         minn = MIN(rr, MIN(gg, bb));
         range = maxx - minn;
         if ((rr > gg) && (rr > bb)) { /* red */
           if (gg > bb) { /* greenish red */
             midrange = rr - gg;
             bb = bb + MIN(YELLOW_DESAT_AMT, range * YELLOW_DESAT_AMT * (midrange/rr) / MAX_DESAT_Y);
           } else { /* blueish red */
             midrange = rr - bb;
             gg = gg + MIN(MAGENTA_DESAT_AMT, range * MAGENTA_DESAT_AMT * (midrange/rr) / MAX_DESAT_M); 
           } /* greenish vs blueish */
           desat_inc = MIN(RED_DESAT_AMT, (midrange * range * RED_DESAT_AMT) / (MAX_DESAT_R * MAX_DESAT_R));
           gg = gg + desat_inc;
           bb = bb + desat_inc;
         } else { /* not red */
           if ((gg > rr) && (gg > bb)) { /* green */
             if (rr > bb) { /* redish green */
               midrange = gg - rr;
               bb = bb + MIN(YELLOW_DESAT_AMT, range * YELLOW_DESAT_AMT * (midrange/gg) / MAX_DESAT_Y);
             } else { /* bluish green */
               midrange = gg - bb;
               rr = rr + MIN(CYAN_DESAT_AMT, range * CYAN_DESAT_AMT * (midrange/gg) / MAX_DESAT_C); 
             } /* redish vs bluish */
             desat_inc = MIN(GRN_DESAT_AMT, (midrange * range * GRN_DESAT_AMT) / (MAX_DESAT_G * MAX_DESAT_G));
             rr = rr + desat_inc;
             bb = bb + desat_inc;
           } else { /* not green */
             if ((bb > rr) && (bb > gg)) { /* blue */
               if (rr > gg) { /* redish blue */
                 midrange = bb - rr;
                 gg = gg + MIN(MAGENTA_DESAT_AMT, range * MAGENTA_DESAT_AMT * (midrange/bb) / MAX_DESAT_M); 
               } else { /* greenish blue */
                 midrange = bb - gg;
                 rr = rr + MIN(CYAN_DESAT_AMT, range * CYAN_DESAT_AMT * (midrange/bb) / MAX_DESAT_C); 
               } /* redish vs greenish blue */
               desat_inc = MIN(BLU_DESAT_AMT, (midrange * range * BLU_DESAT_AMT) / (MAX_DESAT_B * MAX_DESAT_B));
               rr = rr + desat_inc;
               gg = gg + desat_inc;
             } else { /* not blue */
               if ((rr != gg) || (gg != bb) || (rr != bb)) { /* not r=g=b=neutral */
                 if (rr == gg) { /* yellow */
                   bb = bb + MIN(YELLOW_DESAT_AMT,    range * YELLOW_DESAT_AMT  / MAX_DESAT_Y);
                 } else { /* not yellow */
                   if (rr == bb) { /* magenta */
                     gg = gg + MIN(MAGENTA_DESAT_AMT, range * MAGENTA_DESAT_AMT / MAX_DESAT_M); 
                   } else { /* not magenta, is cyan */
                     rr = rr + MIN(CYAN_DESAT_AMT,    range * CYAN_DESAT_AMT    / MAX_DESAT_C); 
                   } /* magenta vs cyan */
                 } /* yellow or not */
               } /* not neutral */
             } /* blue or not */
           } /* green nor not */
         } /* red or not */



          r = (rr >> 2); /* 10bit to 8bit */
          g = (gg >> 2); /* 10bit to 8bit */
          b = (bb >> 2); /* 10bit to 8bit */

          red_frac = rr & 3; /* red fractional low 2 bits */
          grn_frac = gg & 3; /* grn fractional low 2 bits */
          blu_frac = bb & 3; /* blu fractional low 2 bits */

          red_frac_flt = (1.0 * red_frac) / 4.0;
          grn_frac_flt = (1.0 * grn_frac) / 4.0;
          blu_frac_flt = (1.0 * blu_frac) / 4.0;

          r_plus1 = MIN(255, r+1);
          g_plus1 = MIN(255, g+1);
          b_plus1 = MIN(255, b+1);

/* interpolate values using 2lsbs */
#define get_corner(rrr, ggg, bbb, colr) \
          val_indx = rrr * LOOKUP_DIMENSION * LOOKUP_DIMENSION + ggg * LOOKUP_DIMENSION + bbb; \
          colr[0] = lookup[0][val_indx]; \
          colr[1] = lookup[1][val_indx]; \
          colr[2] = lookup[2][val_indx];

          get_corner(r,       g,       b,       col000)
          get_corner(r,       g,       b_plus1, col001)
          get_corner(r,       g_plus1, b,       col010)
          get_corner(r,       g_plus1, b_plus1, col011)
          get_corner(r_plus1, g,       b,       col100)
          get_corner(r_plus1, g,       b_plus1, col101)
          get_corner(r_plus1, g_plus1, b,       col110)
          get_corner(r_plus1, g_plus1, b_plus1, col111)

/* interpolate blue */
          col00[0] = (1.0 - blu_frac_flt) * col000[0] + blu_frac_flt * col001[0];
          col00[1] = (1.0 - blu_frac_flt) * col000[1] + blu_frac_flt * col001[1];
          col00[2] = (1.0 - blu_frac_flt) * col000[2] + blu_frac_flt * col001[2];
          col01[0] = (1.0 - blu_frac_flt) * col010[0] + blu_frac_flt * col011[0];
          col01[1] = (1.0 - blu_frac_flt) * col010[1] + blu_frac_flt * col011[1];
          col01[2] = (1.0 - blu_frac_flt) * col010[2] + blu_frac_flt * col011[2];
          col10[0] = (1.0 - blu_frac_flt) * col100[0] + blu_frac_flt * col101[0];
          col10[1] = (1.0 - blu_frac_flt) * col100[1] + blu_frac_flt * col101[1];
          col10[2] = (1.0 - blu_frac_flt) * col100[2] + blu_frac_flt * col101[2];
          col11[0] = (1.0 - blu_frac_flt) * col110[0] + blu_frac_flt * col111[0];
          col11[1] = (1.0 - blu_frac_flt) * col110[1] + blu_frac_flt * col111[1];
          col11[2] = (1.0 - blu_frac_flt) * col110[2] + blu_frac_flt * col111[2];
/* interpolate green */
          col0[0]  = (1.0 - grn_frac_flt) * col00[0]  + grn_frac_flt * col01[0];
          col0[1]  = (1.0 - grn_frac_flt) * col00[1]  + grn_frac_flt * col01[1];
          col0[2]  = (1.0 - grn_frac_flt) * col00[2]  + grn_frac_flt * col01[2];
          col1[0]  = (1.0 - grn_frac_flt) * col10[0]  + grn_frac_flt * col11[0];
          col1[1]  = (1.0 - grn_frac_flt) * col10[1]  + grn_frac_flt * col11[1];
          col1[2]  = (1.0 - grn_frac_flt) * col10[2]  + grn_frac_flt * col11[2];
/* interpolate red */
          red =      (1.0 - red_frac_flt) * col0[0]   + red_frac_flt * col1[0];
          grn =      (1.0 - red_frac_flt) * col0[1]   + red_frac_flt * col1[1];
          blu =      (1.0 - red_frac_flt) * col0[2]   + red_frac_flt * col1[2];

          pixels_out_exr[y][x].r = red;
          pixels_out_exr[y][x].g = grn;
          pixels_out_exr[y][x].b = blu;
        } /* x */
      } /* y */


    RgbaOutputFile file (outfile, h_reso, v_reso, WRITE_RGB, 1 /* one thread for writing */ );

    file.setFrameBuffer (&pixels_out_exr[0][0], 1, h_reso);
    file.writePixels (v_reso);
}

 printf(" finished writing %s, freeing dpx pixels \n", outfile);
    free(pixels_read);


 } /* frame loop */

} /* main */

