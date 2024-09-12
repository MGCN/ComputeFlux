from paraview.simple import *
import os, sys

inputName= sys.argv[1]
outputName=inputName[0:-6]
outputName=outputName+'.csv'

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
plotOverLine1.Source.Point1 = [0.0, 0.1, 0.0]
plotOverLine1.Source.Point2 = [0.9, 1.0, 0.0]
#plotOverLine1.Source.Point1 = [0.5, 0.0, 0.0]
#plotOverLine1.Source.Point2 = [0.5, 1.0, 0.0]

reader.PointVariables = ['pressure']
SaveData(outputName, proxy=plotOverLine1, Precision=12, WriteTimeSteps=1)
