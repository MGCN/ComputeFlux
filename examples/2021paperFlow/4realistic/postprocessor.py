from paraview.simple import *
import os, sys

inputName=sys.argv[1]
result = inputName.find('_')

pressureHorzLineName=inputName[0:result]+'_H'+inputName[result:len(inputName)]
pressureVertLineName=inputName[0:result]+'_V'+inputName[result:len(inputName)]

pressureHorzLineName=pressureHorzLineName[0:-6]+'.csv'
pressureVertLineName=pressureVertLineName[0:-6]+'.csv'

print(pressureHorzLineName)
print(pressureVertLineName)

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

plotOverLineHorz = PlotOverLine(Input=reader, Source='Line')
plotOverLineHorz.Source.Resolution = 20001
# A line
plotOverLineHorz.Source.Point1 = [  0.0, 500.0, 0.0]
plotOverLineHorz.Source.Point2 = [700.0, 500.0, 0.0]

reader.PointVariables = ['pressure']
SaveData(pressureHorzLineName, proxy=plotOverLineHorz, Precision=12, WriteTimeSteps=1)



plotOverLineVert = PlotOverLine(Input=reader, Source='Line')
plotOverLineVert.Source.Resolution = 20001
# A line
plotOverLineVert.Source.Point1 = [625.0,   0.0, 0.0]
plotOverLineVert.Source.Point2 = [625.0, 600.0, 0.0]

reader.PointVariables = ['pressure']
SaveData(pressureVertLineName, proxy=plotOverLineVert, Precision=12, WriteTimeSteps=1)
