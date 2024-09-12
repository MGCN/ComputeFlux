nx=$1
ref=$2
case=$3

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

if [ "${case}" == "C" ]
then
	perm='1e4'
fi

if [ "${case}" == "B" ]
then
	perm='1e-4'
fi


time mpirun -n ${np} ../../../parrot2-opt -i Diffusion.i nxMoose=${nx} refMoose=${ref} caseMoose=${case} permMoose=${perm}
/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor_${case}.py DiffusionOut_${case}_${nx}_${ref}.e.${np}.0

rm *_0.csv

