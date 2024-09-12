from paraview.simple import *
import os, sys

inputName=sys.argv[1]
case=inputName[-10:-6]
print(case)

if len(sys.argv) > 2:
    outputName = sys.argv[2]
else:
    outputName = 'pres_line_'

pressureLineAName='./'+outputName+case+'.csv'

reader = ExodusIIReader(FileName=inputName)
#Show(ExodusIIReader(FileName=inputName))
#reader = GetActiveSource()
#tsteps = reader.TimestepValues
#print(tsteps)

reader.GenerateObjectIdCellArray = 0
reader.GenerateGlobalElementIdArray = 0
reader.GenerateGlobalNodeIdArray = 0
reader.GlobalVariables = []

#view = GetActiveView()
#view.ViewTime = tsteps[1]
#Render()

plotOverLine1 = PlotOverLine(Input=reader, Source='Line')
plotOverLine1.Source.Resolution = 20001
# A line
plotOverLine1.Source.Point1 = [0.0, 1.1, 0.0]
plotOverLine1.Source.Point2 = [1.0, 1.1, 0.0]

reader.PointVariables = ['pressure']
SaveData(pressureLineAName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)
