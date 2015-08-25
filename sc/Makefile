

	
all: sigma_compare	sigma_compare_PQ
      
sigma_compare : sigma_compare.cpp
	g++ sigma_compare.cpp ../dpxexr/dpx_file_io.cpp -o sigma_compare -lpthread -I /usr/include/OpenEXR -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O1

	
sigma_compare_PQ : sigma_compare_PQ.cpp
	g++ sigma_compare_PQ.cpp ../dpxexr/dpx_file_io.cpp -o sigma_compare_PQ -lpthread -I /usr/include/OpenEXR -lHalf -lIlmImf -lIex -lIlmThread -lz  -m64 -O1
	
clean : 
	rm -v sigma_compare sigma_compare_PQ
  
