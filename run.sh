echo 'running on' $(hostname --short)
cd src/parallel
cp -r make/$(hostname --short) Makefile
make clobber
make
cd ../..

python setup.py

if   [ $(hostname --short) == 'theinfosphere' ] 
then
   #mpirun -np 2 xterm -e gdb ./vhone
   mpirun -np 16 ./vhone
elif [ $(hostname --short) == 'kepler'        ] 
then
   mpirun -np 32 --bind-to core ./vhone
elif [ $(hostname --short) == 'tycho'         ] 
then
   mpirun -np 32 --bind-to-core ./vhone
fi

cd output
./merge_yy
