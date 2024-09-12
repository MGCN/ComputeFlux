nx=$1                  
ny=$(( ${nx} * 225 )) 
ny=$(( ${ny} / 100 ))

np=6

nxString=$nx;

if [ "${nxString}" -lt 100 ]
then
	nxString=0${nxString}
fi

if [ "${nxString}" -lt 1000 ]
then
	nxString=0${nxString}
fi

time mpirun -n ${np} ../../parrot2-opt -i smallFeatures2D.i        nxMoose=${nx} nyMoose=${ny} nxStringMoose=${nxString}
# time mpirun -n ${np} ../../parrot2-opt -i smallFeatures2D_GMG.i    nxMoose=${nx} nyMoose=${ny} nxStringMoose=${nxString}
# time mpirun -n ${np} ../../parrot2-opt -i smallFeatures2D_SSFBBG.i nxMoose=${nx} nyMoose=${ny} nxStringMoose=${nxString}

/Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_${nxString}.e.${np}.0
# /Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_GMG_${nxString}.e.${np}.0    GMG_
# /Applications/ParaView-5.9.0-RC4.app/Contents/bin/pvpython postprocessor.py output_SSFBBG_${nxString}.e.${np}.0 SSFBBG_

rm *_0.csv
