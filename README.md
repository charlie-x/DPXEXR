Tools to help with HDR file manipulation (or rendering) provided by the Gary Demos as listed in the source code. Requires OpenEXR. 

sc/
sigma_compare:  Tool to compare 2 linear (or ordered set) of files that provides a count of pixels at up to 32 stops of range. Computes mean square error between the files ranged to each stop and provides an indication of the width of one-sigma or more of error between the files directly and as a percent of the average value of each stop.  More documentation has been published by Gary through SMPTE.

sigma_compare_PQ: Variations I've made to sigma_compare to provide additional printout of HDR "PQ" 12 bit codevalues based on the floating point (after correction by Amplitude Factor) pixel, sigma and average values.

dpxexr/
Tools to flip files between EXR and DPX32 by Gary.  DPX is a simplier file structure and may be a better mezzanine level uncompressed linear interchange format than 16 bit Integer PQ. 


WARNING: Do not rely on header information of the files processed by these tools. Use at own risk, read the comments and disclaimers in the files. Basically it is just not written completely. Just directly use the floating point values in the files as you intend.
