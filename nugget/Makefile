

	
all: HDR_CPU	
      
HDR_CPU : HDR_CPU.cpp
	g++ HDR_CPU.cpp  -o HDR_CPU -lpthread -I /usr/local/include/OpenEXR -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O3

	
clean : 
	rm -v HDR_CPU 
  
