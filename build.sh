set -x

cd dpxexr
make clean
make

cd ../sc
make clean
make

cd ../nugget
make clean
make

cd ..

pwd 

ls -l dpxexr

ls -l sc

ls -l nugget

