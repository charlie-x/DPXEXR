
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
// This code module reads ACES input files, converts them using the Rendering Nugget and Device Transform into device-specific pixels, and then outputs to either dpx or exr files.


//
// to build:
// linux (tested on Scientific Linux 6.6):
g++ HDR_CPU.cpp -o HDR_CPU_LINUX -lpthread \
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

#define MULTIPLE_THREADS

#ifdef MULTIPLE_THREADS
 #include <pthread.h>

 #define MAXIMUM_NUMBER_OF_THREADS 16

#endif /* MULTIPLE_THREADS */

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

// DPX file reader and writer (only writer is used by this SDK)
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


/*****************************************************************************************************************/
#ifdef MULTIPLE_THREADS
struct rrt_odt_args_struct {
  short rrt_odt_y0, rrt_odt_y1, rrt_odt_h_reso, rrt_odt_v_reso, rrt_odt_odt_type, rrt_odt_thread_num;
  float *rrt_odt_rrt_input, *rrt_odt_odt_result;
}; /* rrt_odt_args_struct */

void
*process_rrt_odt_thread(void *rrt_odt_args_instance)
{
  struct rrt_odt_args_struct *rrt_odt = (struct rrt_odt_args_struct *) rrt_odt_args_instance;
  short y0 = rrt_odt->rrt_odt_y0, y1 = rrt_odt->rrt_odt_y1;
  short thread_num = rrt_odt->rrt_odt_thread_num;
  short h_reso = rrt_odt->rrt_odt_h_reso,  v_reso = rrt_odt->rrt_odt_v_reso;
  short odt_type = rrt_odt-> rrt_odt_odt_type;
  float *rrt_input = rrt_odt->rrt_odt_rrt_input, *odt_result = rrt_odt->rrt_odt_odt_result;

#else /* not MULTIPLE_THREADS */
void
process_rrt_odt(float *rrt_input, float *odt_result, short h_reso, short v_reso, short odt_type)
{

short y0=0;
short y1 = v_reso-1;

#endif /* MULTIPLE_THREADS or not */


short x,y;

float rgbOut[3], aces[3];
int i;


float RGB[3], RGBo[3], RGBoces[3], XYZ[3];

float rOut, gOut, bOut;
float X, Y, Z, R, G, B;
short ans_ready;


// printf(" begin process_rrt_odt for odt_type = %d\n", odt_type);

  for(y=y0; y<y1; y++) {
    for(x=0; x<h_reso; x++) {
      aces[0] = rrt_input[(0*v_reso + y) * h_reso + x];
      aces[1] = rrt_input[(1*v_reso + y) * h_reso + x];
      aces[2] = rrt_input[(2*v_reso + y) * h_reso + x];

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

      odt_result[(0*v_reso + y) * h_reso + x] = rOut;
      odt_result[(1*v_reso + y) * h_reso + x] = gOut;
      odt_result[(2*v_reso + y) * h_reso + x] = bOut;

    } /* x */
  } /* y */


// printf(" end process_rrt_odt\n");

} /* process_rrt_odt */

/***********************************************************************************************************/
int
main(int argc, char **argv)
{
short odt_type=0;
float *odt_result = NULL, *exr_in_pixels = NULL;
short x,y,c,i;
char outfile[300], infile[300];
short num_chars;
int first, last, frame;

short h_reso, v_reso, ans_ready;
short open_exr_vs_dpx_out = 0; /* 0 is open_exr, 1 is dpx */
float s, t, tmp;
int ii;
float expvalue;

#ifdef MULTIPLE_THREADS

  short max_threads = MAXIMUM_NUMBER_OF_THREADS;

  struct rrt_odt_args_struct rrt_odt_args[MAXIMUM_NUMBER_OF_THREADS];
  pthread_t threads[MAXIMUM_NUMBER_OF_THREADS];
  void *stat;

  short y0, y1, ydelta;

#endif /* MULTIPLE_THREADS */

if (argc < 6) {
  printf(" usage: %s infiles, outfiles, first_frame, last_frame, odt_type(1=GD10_Rec709_MDR, 2=GD10_p3_d60_HDR)\n", argv[0]);
  exit(1);
}

 first = atoi(argv[3]);
 last  = atoi(argv[4]);

 odt_type = atoi(argv[5]);

 printf(" processing frames %d to %d with odt_type = %d\n", first, last, odt_type);

#ifdef MULTIPLE_THREADS
  max_threads = MAXIMUM_NUMBER_OF_THREADS;
  printf(" setting maximum number of threads = %d\n", max_threads);
#else /* not MULTIPLE_THREADS */
  max_threads = 1;
  printf(" this version of process_rrt_odt_cpu_multi_threaded is single-threaded, num_threads ignored\n");
#endif /* MULTIPLE_THREADS or not */

   Array2D<Rgba> half_float_pixels;


   for (frame=first; frame <= last; frame++) {


     sprintf(infile, argv[1], frame);
     num_chars = strlen(infile); /* length of outfile string */
     if ((!strcmp(&infile[num_chars-1], "r"))||(!strcmp(&infile[num_chars-1], "R"))) { /* EXR file ending in ".exr" */
      printf(" processing input file %s\n", infile);
     } else { /* not exr */
      printf(" unknown filetype for reading, since extension doesn't end in r, only exr reading supported, infile = %s, aborting\n", outfile);
      exit(1);
     } /* exr or not */

     sprintf(outfile, argv[2], frame);
     num_chars = strlen(outfile); /* length of outfile string */
     if ((!strcmp(&outfile[num_chars-1], "x"))||(!strcmp(&outfile[num_chars-1], "X"))) { /* DPX file ending in ".dpx" */
      printf(" processing output file %s\n", outfile);
      open_exr_vs_dpx_out = 1; /* dpx */
     } else { /* not dpx */
       if ((!strcmp(&outfile[num_chars-1], "r"))||(!strcmp(&outfile[num_chars-1], "R"))) { /* EXR file ending in ".exr" */
        printf(" processing output file %s\n", outfile);
        open_exr_vs_dpx_out = 0; /* exr */
       } else { /* not exr */
        printf(" unknown filetype for writing, since extension doesn't end in x or r, only dpx or exr writing supported, outfile = %s, aborting\n", outfile);
        exit(1);
       } /* exr or not */
     } /* dpx or not */


#ifdef MULTIPLE_THREADS
    RgbaInputFile file (infile, MIN(max_threads , 8 /* arbitrary limit, set as needed */) );
#else /* not MULTIPLE_THREADS */
    RgbaInputFile file (infile, 1 );
#endif /* MULTIPLE_THREADS or not */

    Box2i dw = file.dataWindow();

    h_reso = dw.max.x - dw.min.x + 1;
    v_reso = dw.max.y - dw.min.y + 1;
    half_float_pixels.resizeErase (v_reso, h_reso);

    file.setFrameBuffer (&half_float_pixels[0][0] - dw.min.x - dw.min.y * h_reso, 1, h_reso);
    file.readPixels (dw.min.y, dw.max.y);

    if (exr_in_pixels == NULL) { exr_in_pixels = (float *) malloc(h_reso * v_reso * 12); /* 4-bytes/float * 3-colors */ }

    for(y=0; y< v_reso; y++) {
      for(x=0; x< h_reso; x++) {
        exr_in_pixels[(0*v_reso + y) * h_reso + x] = half_float_pixels[y][x].r; 
        exr_in_pixels[(1*v_reso + y) * h_reso + x] = half_float_pixels[y][x].g; 
        exr_in_pixels[(2*v_reso + y) * h_reso + x] = half_float_pixels[y][x].b; 
      } /* x */
    } /* y */

    if (odt_result == NULL) { odt_result = (float *) malloc(h_reso * v_reso * 12); }

#ifdef MULTIPLE_THREADS
    printf(" perform multi-threaded rrt_odt\n");
    y0 = 0;
    ydelta = MAX(1, v_reso / max_threads);
    y1 = ydelta;
    if (ydelta == 1) { max_threads = MIN( max_threads, v_reso); }

    for (i=0; i< max_threads; i++) {

//      printf(" perform rrt_odt, thread = %d \n", i);

      rrt_odt_args[i].rrt_odt_y0 = y0;
      rrt_odt_args[i].rrt_odt_y1 = y1;
      rrt_odt_args[i].rrt_odt_thread_num = i;
      rrt_odt_args[i].rrt_odt_h_reso  = h_reso;
      rrt_odt_args[i].rrt_odt_v_reso  = v_reso;
      rrt_odt_args[i].rrt_odt_odt_type = odt_type;
      rrt_odt_args[i].rrt_odt_rrt_input = exr_in_pixels;
      rrt_odt_args[i].rrt_odt_odt_result = odt_result;

      pthread_create(&threads[i], NULL, process_rrt_odt_thread, &rrt_odt_args[i]);

      y0 = y1;
      y1 = y1 + ydelta;
      if (i == max_threads - 2) { y1 = v_reso; }

    } /* i */

    for (i=0; i< max_threads; i++) {

      pthread_join(threads[i], &stat);

    } /* i */

   printf(" finished multi-threaded rrt_odt\n");

#else /* not MULTIPLE_THREADS */

  process_rrt_odt(exr_in_pixels, odt_result, h_reso, v_reso, odt_type);

#endif /* MULTIPLE_THREADS or not */

 printf(" begin writing %s\n", outfile);
 if(open_exr_vs_dpx_out == 1) { /* dpx */

   dpx_write_10bit_from_float(outfile, odt_result /*rrt_result*/ /*test_float_ACES_gray_ramp*/, h_reso, v_reso);

 } else { /* exr */

    Array2D<Rgba> pixels_out_exr(v_reso, h_reso);

    for(c=0; c<3; c++) {
      for (y=0; y<v_reso; y++) {
        for (x=0; x<h_reso; x++) {
         tmp = odt_result[(c * v_reso +y) * h_reso +x];
          if (c==0) pixels_out_exr[y][x].r = tmp;
          if (c==1) pixels_out_exr[y][x].g = tmp;
          if (c==2) pixels_out_exr[y][x].b = tmp;
        } /* x */
      } /* y */
    } /* c */

    RgbaOutputFile file (outfile, h_reso, v_reso, WRITE_RGB, 1 /* one thread for writing */ );

    file.setFrameBuffer (&pixels_out_exr[0][0], 1, h_reso);
    file.writePixels (v_reso);

 } /* dpx vs exr */
 printf(" finished writing %s\n", outfile);

 } /* frame loop */

} /* main */

