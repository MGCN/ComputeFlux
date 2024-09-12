nx=$1
ny=$(( ${nx} * 7 )) 
ny=$(( ${ny} / 6 ))
ref=$2

np=6

if [ "${nx}" -lt 10 ]
then
	nx=0${nx}
fi

if [ "${nx}" -lt 100 ]
then
	nx=0${nx}
fi

if [ "${nx}" -lt 1000 ]
then
	nx=0${nx}
fi

if [ "${ref}" -lt 10 ]
then
	ref=0${ref}
fi

time mpirun -n ${np} ../../../parrot2-opt -i Diffusion.i nxMoose=${nx} nyMoose=${ny} refMoose=${ref}
/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py DiffusionOut_${nx}_${ref}.e.${np}.0
rm *_0.csv

#/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython exodus2csv.py DiffusionOut_${nx}_${ref}.e.${np}.0

