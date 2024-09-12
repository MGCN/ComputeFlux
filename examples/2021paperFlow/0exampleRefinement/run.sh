nx=$1
ref=$2

np=1

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

if [ "${refString}" -lt 10 ]
then
	refString=0${refString}
fi

time ../../../parrot2-opt -i createMesh.i nxStringMoose=${nxString} refStringMoose=${refString}

/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython exodus2png.py mesh_${nxString}_${refString}.e mesh_${nxString}_${refString}.jpeg
