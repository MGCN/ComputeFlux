# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------
import os, sys

# print 'Number of arguments:', len(sys.argv), 'arguments.'
# print 'Argument List:', str(sys.argv)

inputName=sys.argv[1]
outputName=inputName[0:14]
outputName = outputName.replace(".e", ".csv")
print(outputName)

folder = os.getcwd()
masterfiles = []
slavefiles = []
name = inputName
for file in os.listdir(folder):
    if file.find(name)>-1:
        masterfiles.append(os.path.join(folder, file))

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

##################
# read master info
##################

# create a new 'ExodusIIReader'
mastere1 = ExodusIIReader(FileName=sorted(masterfiles))

outputNameAndFolder=folder+'/'+outputName;
print(outputNameAndFolder)
writer = CreateWriter(outputNameAndFolder, mastere1, WriteTimeSteps=1)
writer.UpdatePipeline()
