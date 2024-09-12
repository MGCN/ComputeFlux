ref=$1

np=6

refString=$ref;

if [ "${refString}" -lt 10 ]
then
	refString=0${refString}
fi

time mpirun -n ${np} ../../../../parrot2-opt -i createMesh.i refStringMoose=${refString}

/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython exodus2png.py mesh_${refString}.e mesh_${refString}.png
