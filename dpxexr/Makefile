
#
# build dpx/exr tools
#



	
all: EXR_TO_DPX32  DPX32_TO_EXR
      
EXR_TO_DPX32 : EXR_TO_DPX32.cpp
	g++ EXR_TO_DPX32.cpp dpx_file_io.cpp -o EXR_TO_DPX32 -lpthread \
        -I /usr/include/OpenEXR \
         -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O3

DPX32_TO_EXR : DPX32_TO_EXR.cpp
	g++ DPX32_TO_EXR.cpp dpx_file_io.cpp -o DPX32_TO_EXR -lpthread \
 -I /usr/include/OpenEXR \
 -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O3	

	
clean : 
	rm -v EXR_TO_DPX32  DPX32_TO_EXR
  

