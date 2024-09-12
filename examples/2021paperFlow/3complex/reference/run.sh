np=6

time mpirun -n ${np} ../../../../parrot2-opt -i Diffusion_R.i

/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython ../postprocessor.py DiffusionOut_R.e.${np}.0

rm *_0.csv
