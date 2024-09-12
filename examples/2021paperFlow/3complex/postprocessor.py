from paraview.simple import *
import os, sys

inputName= sys.argv[1]
outputName=inputName[0:-6]
outputName=outputName+'.csv'

print(outputName)

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
plotOverLineHorz.Source.Point1 = [0.0, 0.5, 0.0]
plotOverLineHorz.Source.Point2 = [1.0, 0.9, 0.0]

reader.PointVariables = ['pressure']
SaveData(outputName, proxy=plotOverLineHorz, Precision=12, WriteTimeSteps=1)
