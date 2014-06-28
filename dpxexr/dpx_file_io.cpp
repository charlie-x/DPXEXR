
// Author and copyright: Gary Demos 2010,2011,2012,2013,2014
//
//////////////////////////////////////////////////////////////////////////////

/* no warranties expressed nor implied */
/* no representation is herein made as to usefulness nor suitability for any pupose */
/* caution, code may contain bugs and design flaws */
/* use at your own risk */


//
// to build:
// linux:
// g++ dpx_file_io.cpp -o dpx_file_io_LINUX -lpthread -O3 -m64

// note: -O3 will make the code run faster, but some compilers generate bad executables with -O3, -O2, or even -O1.  Most compilers default to -O1

// No attempt is made in this code to write meaningful dpx header information (other than a generic header).
// Programs which rely on dpx header information may thus mis-interpret the generic dpx header which is written by this code.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>

#define DPX_HALF_SUPPORTED
#ifdef DPX_HALF_SUPPORTED

#include <ImfInputFile.h>
#include <ImfRgbaFile.h>
#include <ImfArray.h>

/*#include <ImfBoxAttribute.h>
#include <ImfChannelListAttribute.h>
#include <ImfCompressionAttribute.h>
#include <ImfFloatAttribute.h>
#include <ImfIntAttribute.h>
#include <ImfLineOrderAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfStringAttribute.h>
#include <ImfVecAttribute.h>
#include <iostream>*/

/*using namespace std;
using std::cout;
using std::endl;*/

using namespace Imf;
using namespace Imath;

#endif /* DPX_HALF_SUPPORTED */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define INT_SW(A) (((A >> 24) & 255) | (((A >> 16) & 255) << 8) | (((A >> 8) & 255) << 16) | ((A & 255) << 24));
#define SHORT_SW(A) (((A >> 8) & 255) | ((A & 255) << 8));

#define INTEL_LE

// globals:


// DPX file reader and writer
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
dpx_read (char *inname, float **pixels_read, short *width, short *height, short cineon, short half_flag)
{

FILE *fp;
float *pixels, *pxls_copy;
short x,y;
float red, grn, blu;
short rr, gg, bb;
unsigned int *pixelbuf;
float *pxlbuf;
unsigned short *pxbuf=NULL;
unsigned int *pp;
unsigned int tmp, hdr_offset;
short wide, tall, ww, hh;
int cinhdr[64];
char dpx_header[2048];
short swap_byte_order;
int yvalr, yvalg, yvalb, yval1, yval2, yval3, yv1, yv2, yv3;
short dpx_in_16_or_32 = 0;

  swap_byte_order = 0;

#ifndef DPX_HALF_SUPPORTED
   if (half_flag == 1) { /* convert half-floats to floats */
     printf(" %s half-float reading not supported for dpx files in this version, aborting\n", inname);
     exit(1);
   } /* half_flag */
#endif /* not DPX_HALF_SUPPORTED */

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
      dpx_in_16_or_32 = 2; /* 2 is 32bit */
#ifdef VERBOSE_DECODE
      printf(" dpx packing is 32-bit float \n ");
#endif /* VERBOSE_DECODE */
    } else {
      if (dpx_header[803] == 16) {
#ifdef VERBOSE_DECODE
        printf(" dpx packing is 16-bit \n ");
#endif /* VERBOSE_DECODE */
        dpx_in_16_or_32 = 1; /* 1 is 16bit */
      } else { /* not 16bits */
        if (dpx_header[803] == 12) {
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

  if (dpx_in_16_or_32 == 2) { /* floats */
    pxlbuf = (float *) malloc(wide*tall*12); /* rgb float */
  } else { /* not float */
    if (dpx_in_16_or_32 == 1) { /* 16-bits */
      pxbuf = (unsigned short *) malloc(wide*tall*6); /* 2 bytes per color rgb */
    } else { /* not 16bit, must be 10bit */
      pixelbuf = (unsigned int *) malloc(wide*tall*4); /* 4 bytes per pixel in 10-bit dpx packing */
    } /* 16bits or not */
  } /* dpx_in_16_or_32 == 2 or not */



  pixels = (float *) malloc(wide * tall * 12); /* 4-bytes/float * 3-colors */
  *pixels_read = &pixels[0];

#ifdef VERBOSE_PROCESSING 
 printf(" reading dpx or cineon file pixels \n ");
#endif /* VERBOSE_PROCESSING */

  if (dpx_in_16_or_32 == 2) { /* floats */
    fread( pxlbuf, 1, wide * tall * sizeof(float) * 3, fp);
    if (ferror(fp)) {
      printf(" error reading data for file %s in dpx_read float, aborting\n", y, inname);
      exit(1);
    }
  } else { /* not float */
    if (dpx_in_16_or_32 == 1) { /* ushorts */
      fread( pxbuf, 1, wide * tall * sizeof(unsigned short) * 3, fp);
      if (ferror(fp)) {
        printf(" error reading data for file %s in dpx_read float, aborting\n", y, inname);
        exit(1);
      }
    } else { /* not ushorts */
/* must be 10bit */
      for( y = 0; y < tall; y++) {
        pp = &pixelbuf[0] + y * wide;
        fread( pp, 1, wide * sizeof(unsigned int), fp);
        if (ferror(fp)) {
          printf(" error reading data for scanline %d for file %s in dpx_read (10bit dpx file), aborting\n", y, inname);
          exit(1);
        }
      } /* y */
    } /* ushorts or not */
  } /* float or not */

 fclose(fp);

  if (dpx_in_16_or_32 == 2) { /* floats */

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

  } else { /* not float */

   if (dpx_in_16_or_32 == 1) { /* ushort */
/* caution, the 32-bit packing flag in the header can be misleading for ushorts, since the dpx2.0 doc says little endian, whereas this is nonsense.  ushorts will always go in rgb normal order, and have no 32-bit attribute */
    for( y = 0; y < tall; y++) {
      yvalr = (0*tall + y)*wide;
      yvalg = (1*tall + y)*wide;
      yvalb = (2*tall + y)*wide;
      for( x = 0; x < wide ; x++) {

        if (half_flag == 1) { /* convert half-floats to floats */
#ifdef DPX_HALF_SUPPORTED
unsigned short red_tmp, grn_tmp, blu_tmp;
          if (swap_byte_order == 1) {
unsigned short ttmp;
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3]);
            red_tmp = SHORT_SW(ttmp);
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3 + 1]);
            grn_tmp = SHORT_SW(ttmp);
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3 + 2]);
            blu_tmp = SHORT_SW(ttmp);
          } else { /* dont swap byte order */
            red_tmp=pxbuf[(y * wide + x) * 3];
            grn_tmp=pxbuf[(y * wide + x) * 3 + 1];
            blu_tmp=pxbuf[(y * wide + x) * 3 + 2];
          }/* swap_byte_order or not */
half hlf_tmp;
          hlf_tmp = *((half *) &red_tmp);
          red = hlf_tmp; /* half to float */
          hlf_tmp = *((half *) &grn_tmp);
          grn = hlf_tmp; /* half to float */
          hlf_tmp = *((half *) &blu_tmp);
          blu = hlf_tmp; /* half to float */
#endif /* DPX_HALF_SUPPORTED or not */
        } else { /* not half_flag, ushorts in the range 0.0 to 1.0 */
          if (swap_byte_order == 1) {
unsigned short ttmp;
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3]);
            ttmp = SHORT_SW(ttmp);
            red = ttmp / 65535.0;
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3 + 1]);
            ttmp = SHORT_SW(ttmp);
            grn = ttmp / 65535.0;
            ttmp=*((unsigned short *) &pxbuf[(y * wide + x) * 3 + 2]);
            ttmp = SHORT_SW(ttmp);
            blu = ttmp / 65535.0;
          } else { /* dont swap byte order */
            red=pxbuf[(y * wide + x) * 3]     / 65535.0;
            grn=pxbuf[(y * wide + x) * 3 + 1] / 65535.0;
            blu=pxbuf[(y * wide + x) * 3 + 2] / 65535.0;
          }/* swap_byte_order or not */
        } /* half_flag or not */

        pixels[yvalr + x] = red;
        pixels[yvalg + x] = grn;
        pixels[yvalb + x] = blu;

      } /* x */
    } /* y */

   } else { /* not ushort */
/* 10bit */
     for( y = 0; y < tall; y++) {
      yvalr = (0*tall + y)*wide;
      yvalg = (1*tall + y)*wide;
      yvalb = (2*tall + y)*wide;
      for( x = 0; x < wide ; x++) {

        if (swap_byte_order == 1) {
          tmp=INT_SW(pixelbuf[y * wide + x]);
        } else { /* dont swap byte order */
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
  } /* ushort or not */
 } /* float or not */


  if (dpx_in_16_or_32 == 2) { /* floats */
    free(pxlbuf);
   } else { /* not float */
    if (dpx_in_16_or_32 == 1) { /* ushorts */
     free(pxbuf);
    } else { /* 10bit */
     free(pixelbuf);
    } /* ushorts or not */
  } /* float or not */


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
 dpx_header[801] = 2; /* transfer characteristic, 6 = video_gamma, 1 = printing density, 2 = lionear */
 dpx_header[802] = 4; /* colorimetric, 6 = video_gamma, 1 = printing density, 4 = unspecified color */
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
 dpx_header[801] = 2; /* transfer characteristic, 6 = video_gamma, 1 = printing density, 2 = linear */
 dpx_header[802] = 4; /* colorimetric, 6 = video_gamma, 1 = printing density, 4 = unspecified color */
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
/* write non-standard half-float pixels in dpx file */
void
dpx_write_half(char *outname, float *pixel_result, short width, short height)
{

#ifdef DPX_HALF_SUPPORTED

FILE *fp_out;

unsigned int tmp, size;
short x, y, i;
float red, grn, blu;

unsigned short *pixelbuf;
unsigned char dpx_header[HEADER_SIZE];
int yvalr, yvalg, yvalb, yval0;


     if ((fp_out = fopen(outname,"wb" /* rb for reading, wb for writing */)) == NULL) {
        printf(" Cannot open dpx output file %s \n", outname);
        exit(1);
      } /* fopen */

  for (i=0; i< HEADER_SIZE; i++) { /* clear header array */
    dpx_header[i]=0;
  } /* i */


  pixelbuf = (unsigned short *) malloc(width*height*6); /* rgb half-floats */

half half_red, half_grn, half_blu;
 for( y = 0; y < height; y++) {
  yval0 = y * width * 3;
  yvalr = (0*height + y)*width;
  yvalg = (1*height + y)*width;
  yvalb = (2*height + y)*width;
  for( x = 0; x < width ; x++) {

      red = pixel_result[yvalr + x];
      grn = pixel_result[yvalg + x];
      blu = pixel_result[yvalb + x];

      half_red = red; /* float to half */
      half_grn = grn; /* float to half */
      half_blu = blu; /* float to half */
#ifdef INTEL_LE
unsigned short ttmp;
      ttmp = *((unsigned short *) &half_red);
      *((unsigned short *) &half_red) = SHORT_SW(ttmp)
      ttmp = *((unsigned short *) &half_grn);
      *((unsigned short *) &half_grn) = SHORT_SW(ttmp)
      ttmp = *((unsigned short *) &half_blu);
      *((unsigned short *) &half_blu) = SHORT_SW(ttmp)
#endif /* INTEL_LE */

      pixelbuf[yval0 + 3 * x    ] = *((unsigned short *) &half_red);
      pixelbuf[yval0 + 3 * x + 1] = *((unsigned short *) &half_grn);
      pixelbuf[yval0 + 3 * x + 2] = *((unsigned short *) &half_blu);

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

  size = width * height * 6; /* 6 bytes per pixel, (in bytes???) is this right for the size field???, or should this be in pixels??? */

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
 dpx_header[801] = 2; /* transfer characteristic, 6 = video_gamma, 1 = printing density, 2 = linear */
 dpx_header[802] = 4; /* colorimetric, 6 = video_gamma, 1 = printing density, 4 = unspecified color */
 dpx_header[803] = 16; /* bits per element */
 *((short *) &dpx_header[804]) = 0; /* packing, no packing, note: this is confusing, since there is no 32-bit packing with respect to 16-bit rgb values, the dpx2.0 doc says little-endian, but that is nonsense */
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
 dpx_header[801] = 2; /* transfer characteristic, 6 = video_gamma, 1 = printing density, 2 = linear */
 dpx_header[802] = 4; /* colorimetric, 6 = video_gamma, 1 = printing density, 4 = unspecified color */
 dpx_header[803] = 16; /* bits per element */
 *((short *) &dpx_header[804]) = 0; /* packing, no packing, note: this is confusing, since there is no 32-bit packing with respect to 16-bit rgb values, the dpx2.0 doc says little-endian, but that is nonsense */
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

    fwrite( pixelbuf, 1, height*width*6, fp_out);

    if (ferror(fp_out)) {
      printf(" error writing data to file %s in dpx_write, aborting\n", outname);
      exit(1);
    }

#ifdef VERBOSE_DECODE
    printf(" closing dpxhlf half-float (nonstandard) dpx output file \n");
#endif /* VERBOSE_DECODE */

 if(fclose(fp_out)) {
   printf(" error closing file %s in dpx_write_half, aborting\n", outname);
   exit(1);
 }

    free(pixelbuf);
#else /* not DPX_HALF_SUPPORTED */
printf(" %s dpx_half writing not supported in this version, aborting\n", outname);
exit(1);
#endif /* DPX_HALF_SUPPORTED or not */


} /* dpx_write_half */

/***********************************************************************************************************/

#if 0
/***********************************************************************************************************/
int
main(int argc, char **argv)
{
float *pixels = NULL;
short x,y,c,i;
char outfile[300], infile[300];
short num_chars;
int first, last, frame;
short cineon=0;
short h_reso, v_reso;
float s, t, tmp;
int ii;

if (argc < 5) {
  printf(" usage: %s infiles, outfiles, first_frame, last_frame\n", argv[0]);
  exit(1);
}

 first = atoi(argv[3]);
 last  = atoi(argv[4]);

 printf(" processing frames %d to %d\n", first, last);


   for (frame=first; frame <= last; frame++) {


     sprintf(infile, argv[1], frame);
     num_chars = strlen(infile); /* length of outfile string */
     if ((!strcmp(&infile[num_chars-1], "x"))||(!strcmp(&infile[num_chars-1], "X"))) { /* DPX file ending in ".dpx" */
      printf(" processing input file %s\n", infile);
     } else { /* not exr */
      printf(" unknown filetype for reading, since extension doesn't end in x, only dpx reading supported, infile = %s, aborting\n", outfile);
      exit(1);
     } /* exr or not */

     sprintf(outfile, argv[2], frame);
     num_chars = strlen(outfile); /* length of outfile string */
     if ((!strcmp(&outfile[num_chars-1], "x"))||(!strcmp(&outfile[num_chars-1], "X"))) { /* DPX file ending in ".dpx" */
      printf(" processing output file %s\n", outfile);
     } else { /* not dpx */
printf(" output file %s not supported, aborting\n", outfile);
exit(1);
     } /* dpx or not */

    dpx_read (infile, &pixels, &h_reso, &v_reso, cineon /* make cineon nonzero for cineon file reading */);

    dpx_write_10bit_from_float(outfile, pixels, h_reso, v_reso);

 printf(" finished writing %s\n", outfile);

   free(pixels);

 } /* frame loop */

} /* main */
#endif /* 0 */

