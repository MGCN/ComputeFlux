nx=$1
ref=$2
dir=$3

np=6

refString=$ref;
nxString=$nx;

if [ "${nxString}" -lt 10 ]
then
	nxString=0${nxString}
fi

if [ "${nxString}" -lt 100 ]
then
	nxString=0${nxString}
fi

if [ "${nxString}" -lt 1000 ]
then
	nxString=0${nxString}
fi


if [ "${refString}" -lt 10 ]
then
	refString=0${refString}
fi

time mpirun -n ${np} ../../../parrot2-opt -i Diffusion_${dir}.i nxStringMoose=${nxString} refStringMoose=${refString}

/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py DiffusionOut_${dir}_${nxString}_${refString}.e.${np}.0

rm *_0.csv
