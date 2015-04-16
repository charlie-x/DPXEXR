

	
all: HDR_CPU HDR_CPU_INVERSE
      
HDR_CPU : HDR_CPU.cpp hdr_macros.h
	g++ HDR_CPU.cpp  -o HDR_CPU -lpthread -I /usr/local/include/OpenEXR -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O3
	
HDR_CPU_INVERSE : HDR_CPU_INVERSE.cpp hdr_macros.h
	g++ HDR_CPU_INVERSE.cpp  -o HDR_CPU_INVERSE -lpthread -I /usr/local/include/OpenEXR -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O3
	
clean : 
	rm -v HDR_CPU HDR_CPU_INVERSE
  
