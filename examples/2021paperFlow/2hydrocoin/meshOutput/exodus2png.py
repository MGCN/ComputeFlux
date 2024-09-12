# trace generated using paraview version 5.9.0-RC4

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

inputName= sys.argv[1]
outputName=sys.argv[2]


# create a new 'ExodusIIReader'
mesh_01e = ExodusIIReader(registrationName='mesh_01.e', FileName=inputName)
mesh_01e.GenerateObjectIdCellArray = 1
mesh_01e.GenerateGlobalElementIdArray = 1
mesh_01e.ElementVariables = []
mesh_01e.FaceVariables = []
mesh_01e.EdgeVariables = []
mesh_01e.SideSetResultArrayStatus = []
mesh_01e.NodeSetResultArrayStatus = []
mesh_01e.FaceSetResultArrayStatus = []
mesh_01e.EdgeSetResultArrayStatus = []
mesh_01e.GenerateGlobalNodeIdArray = 1
mesh_01e.ElementSetResultArrayStatus = []
mesh_01e.PointVariables = []
mesh_01e.GlobalVariables = []
mesh_01e.ApplyDisplacements = 1
mesh_01e.DisplacementMagnitude = 1.0
mesh_01e.EdgeBlocks = []
mesh_01e.NodeSetArrayStatus = []
mesh_01e.SideSetArrayStatus = []
mesh_01e.FaceSetArrayStatus = []
mesh_01e.EdgeSetArrayStatus = []
mesh_01e.ElementSetArrayStatus = []
mesh_01e.NodeMapArrayStatus = []
mesh_01e.EdgeMapArrayStatus = []
mesh_01e.FaceMapArrayStatus = []
mesh_01e.ElementMapArrayStatus = []
mesh_01e.ElementBlocks = []
mesh_01e.FaceBlocks = []
mesh_01e.HasModeShapes = 0
mesh_01e.ModeShape = 1
mesh_01e.AnimateVibrations = 1
mesh_01e.IgnoreFileTime = 0
mesh_01e.GenerateFileIdArray = 0
mesh_01e.UseLegacyBlockNamesWithElementTypes = 0

# Properties modified on mesh_01e
mesh_01e.ElementBlocks = ['Unnamed block ID: 1']
mesh_01e.FilePrefix = ''
mesh_01e.FilePattern = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
mesh_01eDisplay = Show(mesh_01e, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
mesh_01eDisplay.Selection = None
mesh_01eDisplay.Representation = 'Surface'
mesh_01eDisplay.ColorArrayName = [None, '']
mesh_01eDisplay.LookupTable = None
mesh_01eDisplay.MapScalars = 1
mesh_01eDisplay.MultiComponentsMapping = 0
mesh_01eDisplay.InterpolateScalarsBeforeMapping = 1
mesh_01eDisplay.Opacity = 1.0
mesh_01eDisplay.PointSize = 2.0
mesh_01eDisplay.LineWidth = 1.0
mesh_01eDisplay.RenderLinesAsTubes = 0
mesh_01eDisplay.RenderPointsAsSpheres = 0
mesh_01eDisplay.Interpolation = 'Gouraud'
mesh_01eDisplay.Specular = 0.0
mesh_01eDisplay.SpecularColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.SpecularPower = 100.0
mesh_01eDisplay.Luminosity = 0.0
mesh_01eDisplay.Ambient = 0.0
mesh_01eDisplay.Diffuse = 1.0
mesh_01eDisplay.Roughness = 0.3
mesh_01eDisplay.Metallic = 0.0
mesh_01eDisplay.EdgeTint = [1.0, 1.0, 1.0]
mesh_01eDisplay.SelectTCoordArray = 'None'
mesh_01eDisplay.SelectNormalArray = 'None'
mesh_01eDisplay.SelectTangentArray = 'None'
mesh_01eDisplay.Texture = None
mesh_01eDisplay.RepeatTextures = 1
mesh_01eDisplay.InterpolateTextures = 0
mesh_01eDisplay.SeamlessU = 0
mesh_01eDisplay.SeamlessV = 0
mesh_01eDisplay.UseMipmapTextures = 0
mesh_01eDisplay.BaseColorTexture = None
mesh_01eDisplay.NormalTexture = None
mesh_01eDisplay.NormalScale = 1.0
mesh_01eDisplay.MaterialTexture = None
mesh_01eDisplay.OcclusionStrength = 1.0
mesh_01eDisplay.EmissiveTexture = None
mesh_01eDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
mesh_01eDisplay.FlipTextures = 0
mesh_01eDisplay.BackfaceRepresentation = 'Follow Frontface'
mesh_01eDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.BackfaceOpacity = 1.0
mesh_01eDisplay.Position = [0.0, 0.0, 0.0]
mesh_01eDisplay.Scale = [1.0, 1.0, 1.0]
mesh_01eDisplay.Orientation = [0.0, 0.0, 0.0]
mesh_01eDisplay.Origin = [0.0, 0.0, 0.0]
mesh_01eDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
mesh_01eDisplay.Pickable = 1
mesh_01eDisplay.Triangulate = 0
mesh_01eDisplay.UseShaderReplacements = 0
mesh_01eDisplay.ShaderReplacements = ''
mesh_01eDisplay.NonlinearSubdivisionLevel = 1
mesh_01eDisplay.UseDataPartitions = 0
mesh_01eDisplay.OSPRayUseScaleArray = 'All Approximate'
mesh_01eDisplay.OSPRayScaleArray = 'GlobalNodeId'
mesh_01eDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
mesh_01eDisplay.OSPRayMaterial = 'None'
mesh_01eDisplay.Orient = 0
mesh_01eDisplay.OrientationMode = 'Direction'
mesh_01eDisplay.SelectOrientationVectors = 'None'
mesh_01eDisplay.Scaling = 0
mesh_01eDisplay.ScaleMode = 'No Data Scaling Off'
mesh_01eDisplay.ScaleFactor = 160.0
mesh_01eDisplay.SelectScaleArray = 'GlobalNodeId'
mesh_01eDisplay.GlyphType = 'Arrow'
mesh_01eDisplay.UseGlyphTable = 0
mesh_01eDisplay.GlyphTableIndexArray = 'GlobalNodeId'
mesh_01eDisplay.UseCompositeGlyphTable = 0
mesh_01eDisplay.UseGlyphCullingAndLOD = 0
mesh_01eDisplay.LODValues = []
mesh_01eDisplay.ColorByLODIndex = 0
mesh_01eDisplay.GaussianRadius = 8.0
mesh_01eDisplay.ShaderPreset = 'Sphere'
mesh_01eDisplay.CustomTriangleScale = 3
mesh_01eDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
mesh_01eDisplay.Emissive = 0
mesh_01eDisplay.ScaleByArray = 0
mesh_01eDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
mesh_01eDisplay.ScaleArrayComponent = ''
mesh_01eDisplay.UseScaleFunction = 1
mesh_01eDisplay.ScaleTransferFunction = 'PiecewiseFunction'
mesh_01eDisplay.OpacityByArray = 0
mesh_01eDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
mesh_01eDisplay.OpacityArrayComponent = ''
mesh_01eDisplay.OpacityTransferFunction = 'PiecewiseFunction'
mesh_01eDisplay.DataAxesGrid = 'GridAxesRepresentation'
mesh_01eDisplay.SelectionCellLabelBold = 0
mesh_01eDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
mesh_01eDisplay.SelectionCellLabelFontFamily = 'Arial'
mesh_01eDisplay.SelectionCellLabelFontFile = ''
mesh_01eDisplay.SelectionCellLabelFontSize = 18
mesh_01eDisplay.SelectionCellLabelItalic = 0
mesh_01eDisplay.SelectionCellLabelJustification = 'Left'
mesh_01eDisplay.SelectionCellLabelOpacity = 1.0
mesh_01eDisplay.SelectionCellLabelShadow = 0
mesh_01eDisplay.SelectionPointLabelBold = 0
mesh_01eDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
mesh_01eDisplay.SelectionPointLabelFontFamily = 'Arial'
mesh_01eDisplay.SelectionPointLabelFontFile = ''
mesh_01eDisplay.SelectionPointLabelFontSize = 18
mesh_01eDisplay.SelectionPointLabelItalic = 0
mesh_01eDisplay.SelectionPointLabelJustification = 'Left'
mesh_01eDisplay.SelectionPointLabelOpacity = 1.0
mesh_01eDisplay.SelectionPointLabelShadow = 0
mesh_01eDisplay.PolarAxes = 'PolarAxesRepresentation'
mesh_01eDisplay.ScalarOpacityFunction = None
mesh_01eDisplay.ScalarOpacityUnitDistance = 343.34848133744094
mesh_01eDisplay.UseSeparateOpacityArray = 0
mesh_01eDisplay.OpacityArrayName = ['POINTS', 'GlobalNodeId']
mesh_01eDisplay.OpacityComponent = ''
mesh_01eDisplay.ExtractedBlockIndex = 2
mesh_01eDisplay.SelectMapper = 'Projected tetra'
mesh_01eDisplay.SamplingDimensions = [128, 128, 128]
mesh_01eDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
mesh_01eDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
mesh_01eDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
mesh_01eDisplay.GlyphType.TipResolution = 6
mesh_01eDisplay.GlyphType.TipRadius = 0.1
mesh_01eDisplay.GlyphType.TipLength = 0.35
mesh_01eDisplay.GlyphType.ShaftResolution = 6
mesh_01eDisplay.GlyphType.ShaftRadius = 0.03
mesh_01eDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
mesh_01eDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 237.0, 1.0, 0.5, 0.0]
mesh_01eDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
mesh_01eDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 237.0, 1.0, 0.5, 0.0]
mesh_01eDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
mesh_01eDisplay.DataAxesGrid.XTitle = 'X Axis'
mesh_01eDisplay.DataAxesGrid.YTitle = 'Y Axis'
mesh_01eDisplay.DataAxesGrid.ZTitle = 'Z Axis'
mesh_01eDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.XTitleFontFile = ''
mesh_01eDisplay.DataAxesGrid.XTitleBold = 0
mesh_01eDisplay.DataAxesGrid.XTitleItalic = 0
mesh_01eDisplay.DataAxesGrid.XTitleFontSize = 12
mesh_01eDisplay.DataAxesGrid.XTitleShadow = 0
mesh_01eDisplay.DataAxesGrid.XTitleOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.YTitleFontFile = ''
mesh_01eDisplay.DataAxesGrid.YTitleBold = 0
mesh_01eDisplay.DataAxesGrid.YTitleItalic = 0
mesh_01eDisplay.DataAxesGrid.YTitleFontSize = 12
mesh_01eDisplay.DataAxesGrid.YTitleShadow = 0
mesh_01eDisplay.DataAxesGrid.YTitleOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.ZTitleFontFile = ''
mesh_01eDisplay.DataAxesGrid.ZTitleBold = 0
mesh_01eDisplay.DataAxesGrid.ZTitleItalic = 0
mesh_01eDisplay.DataAxesGrid.ZTitleFontSize = 12
mesh_01eDisplay.DataAxesGrid.ZTitleShadow = 0
mesh_01eDisplay.DataAxesGrid.ZTitleOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.FacesToRender = 63
mesh_01eDisplay.DataAxesGrid.CullBackface = 0
mesh_01eDisplay.DataAxesGrid.CullFrontface = 1
mesh_01eDisplay.DataAxesGrid.ShowGrid = 0
mesh_01eDisplay.DataAxesGrid.ShowEdges = 1
mesh_01eDisplay.DataAxesGrid.ShowTicks = 1
mesh_01eDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
mesh_01eDisplay.DataAxesGrid.AxesToLabel = 63
mesh_01eDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.XLabelFontFile = ''
mesh_01eDisplay.DataAxesGrid.XLabelBold = 0
mesh_01eDisplay.DataAxesGrid.XLabelItalic = 0
mesh_01eDisplay.DataAxesGrid.XLabelFontSize = 12
mesh_01eDisplay.DataAxesGrid.XLabelShadow = 0
mesh_01eDisplay.DataAxesGrid.XLabelOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.YLabelFontFile = ''
mesh_01eDisplay.DataAxesGrid.YLabelBold = 0
mesh_01eDisplay.DataAxesGrid.YLabelItalic = 0
mesh_01eDisplay.DataAxesGrid.YLabelFontSize = 12
mesh_01eDisplay.DataAxesGrid.YLabelShadow = 0
mesh_01eDisplay.DataAxesGrid.YLabelOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
mesh_01eDisplay.DataAxesGrid.ZLabelFontFile = ''
mesh_01eDisplay.DataAxesGrid.ZLabelBold = 0
mesh_01eDisplay.DataAxesGrid.ZLabelItalic = 0
mesh_01eDisplay.DataAxesGrid.ZLabelFontSize = 12
mesh_01eDisplay.DataAxesGrid.ZLabelShadow = 0
mesh_01eDisplay.DataAxesGrid.ZLabelOpacity = 1.0
mesh_01eDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
mesh_01eDisplay.DataAxesGrid.XAxisPrecision = 2
mesh_01eDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
mesh_01eDisplay.DataAxesGrid.XAxisLabels = []
mesh_01eDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
mesh_01eDisplay.DataAxesGrid.YAxisPrecision = 2
mesh_01eDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
mesh_01eDisplay.DataAxesGrid.YAxisLabels = []
mesh_01eDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
mesh_01eDisplay.DataAxesGrid.ZAxisPrecision = 2
mesh_01eDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
mesh_01eDisplay.DataAxesGrid.ZAxisLabels = []
mesh_01eDisplay.DataAxesGrid.UseCustomBounds = 0
mesh_01eDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
mesh_01eDisplay.PolarAxes.Visibility = 0
mesh_01eDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
mesh_01eDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
mesh_01eDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
mesh_01eDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
mesh_01eDisplay.PolarAxes.EnableCustomRange = 0
mesh_01eDisplay.PolarAxes.CustomRange = [0.0, 1.0]
mesh_01eDisplay.PolarAxes.PolarAxisVisibility = 1
mesh_01eDisplay.PolarAxes.RadialAxesVisibility = 1
mesh_01eDisplay.PolarAxes.DrawRadialGridlines = 1
mesh_01eDisplay.PolarAxes.PolarArcsVisibility = 1
mesh_01eDisplay.PolarAxes.DrawPolarArcsGridlines = 1
mesh_01eDisplay.PolarAxes.NumberOfRadialAxes = 0
mesh_01eDisplay.PolarAxes.AutoSubdividePolarAxis = 1
mesh_01eDisplay.PolarAxes.NumberOfPolarAxis = 0
mesh_01eDisplay.PolarAxes.MinimumRadius = 0.0
mesh_01eDisplay.PolarAxes.MinimumAngle = 0.0
mesh_01eDisplay.PolarAxes.MaximumAngle = 90.0
mesh_01eDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
mesh_01eDisplay.PolarAxes.Ratio = 1.0
mesh_01eDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
mesh_01eDisplay.PolarAxes.PolarAxisTitleVisibility = 1
mesh_01eDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
mesh_01eDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
mesh_01eDisplay.PolarAxes.PolarLabelVisibility = 1
mesh_01eDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
mesh_01eDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
mesh_01eDisplay.PolarAxes.RadialLabelVisibility = 1
mesh_01eDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
mesh_01eDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
mesh_01eDisplay.PolarAxes.RadialUnitsVisibility = 1
mesh_01eDisplay.PolarAxes.ScreenSize = 10.0
mesh_01eDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
mesh_01eDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
mesh_01eDisplay.PolarAxes.PolarAxisTitleFontFile = ''
mesh_01eDisplay.PolarAxes.PolarAxisTitleBold = 0
mesh_01eDisplay.PolarAxes.PolarAxisTitleItalic = 0
mesh_01eDisplay.PolarAxes.PolarAxisTitleShadow = 0
mesh_01eDisplay.PolarAxes.PolarAxisTitleFontSize = 12
mesh_01eDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
mesh_01eDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
mesh_01eDisplay.PolarAxes.PolarAxisLabelFontFile = ''
mesh_01eDisplay.PolarAxes.PolarAxisLabelBold = 0
mesh_01eDisplay.PolarAxes.PolarAxisLabelItalic = 0
mesh_01eDisplay.PolarAxes.PolarAxisLabelShadow = 0
mesh_01eDisplay.PolarAxes.PolarAxisLabelFontSize = 12
mesh_01eDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
mesh_01eDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
mesh_01eDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
mesh_01eDisplay.PolarAxes.LastRadialAxisTextBold = 0
mesh_01eDisplay.PolarAxes.LastRadialAxisTextItalic = 0
mesh_01eDisplay.PolarAxes.LastRadialAxisTextShadow = 0
mesh_01eDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
mesh_01eDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
mesh_01eDisplay.PolarAxes.EnableDistanceLOD = 1
mesh_01eDisplay.PolarAxes.DistanceLODThreshold = 0.7
mesh_01eDisplay.PolarAxes.EnableViewAngleLOD = 1
mesh_01eDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
mesh_01eDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
mesh_01eDisplay.PolarAxes.PolarTicksVisibility = 1
mesh_01eDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
mesh_01eDisplay.PolarAxes.TickLocation = 'Both'
mesh_01eDisplay.PolarAxes.AxisTickVisibility = 1
mesh_01eDisplay.PolarAxes.AxisMinorTickVisibility = 0
mesh_01eDisplay.PolarAxes.ArcTickVisibility = 1
mesh_01eDisplay.PolarAxes.ArcMinorTickVisibility = 0
mesh_01eDisplay.PolarAxes.DeltaAngleMajor = 10.0
mesh_01eDisplay.PolarAxes.DeltaAngleMinor = 5.0
mesh_01eDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
mesh_01eDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
mesh_01eDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
mesh_01eDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
mesh_01eDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
mesh_01eDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
mesh_01eDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
mesh_01eDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
mesh_01eDisplay.PolarAxes.ArcMajorTickSize = 0.0
mesh_01eDisplay.PolarAxes.ArcTickRatioSize = 0.3
mesh_01eDisplay.PolarAxes.ArcMajorTickThickness = 1.0
mesh_01eDisplay.PolarAxes.ArcTickRatioThickness = 0.5
mesh_01eDisplay.PolarAxes.Use2DMode = 0
mesh_01eDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [800.0, -425.0, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh_01eDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
mesh_01eDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# change representation type
mesh_01eDisplay.SetRepresentationType('Surface With Edges')

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1982, 1256)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [800.0, -425.0, 10000.0]
renderView1.CameraFocalPoint = [800.0, -425.0, 0.0]
renderView1.CameraParallelScale = 985.2030247619016

# save screenshot
SaveScreenshot(outputName, renderView1, ImageResolution=[1600, 1200],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1982, 1256)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [800.0, -425.0, 10000.0]
renderView1.CameraFocalPoint = [800.0, -425.0, 0.0]
renderView1.CameraParallelScale = 985.2030247619016

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).