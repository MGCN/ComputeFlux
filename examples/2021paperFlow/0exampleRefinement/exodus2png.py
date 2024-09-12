# trace generated using paraview version 5.9.0-RC4

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

inputName= sys.argv[1]
outputName=sys.argv[2]

# create a new 'ExodusIIReader'
mesh_010_08e = ExodusIIReader(registrationName='mesh_010_08.e', FileName=inputName)
mesh_010_08e.GenerateObjectIdCellArray = 1
mesh_010_08e.GenerateGlobalElementIdArray = 1
mesh_010_08e.ElementVariables = []
mesh_010_08e.FaceVariables = []
mesh_010_08e.EdgeVariables = []
mesh_010_08e.SideSetResultArrayStatus = []
mesh_010_08e.NodeSetResultArrayStatus = []
mesh_010_08e.FaceSetResultArrayStatus = []
mesh_010_08e.EdgeSetResultArrayStatus = []
mesh_010_08e.GenerateGlobalNodeIdArray = 1
mesh_010_08e.ElementSetResultArrayStatus = []
mesh_010_08e.PointVariables = []
mesh_010_08e.GlobalVariables = []
mesh_010_08e.ApplyDisplacements = 1
mesh_010_08e.DisplacementMagnitude = 1.0
mesh_010_08e.EdgeBlocks = []
mesh_010_08e.NodeSetArrayStatus = []
mesh_010_08e.SideSetArrayStatus = []
mesh_010_08e.FaceSetArrayStatus = []
mesh_010_08e.EdgeSetArrayStatus = []
mesh_010_08e.ElementSetArrayStatus = []
mesh_010_08e.NodeMapArrayStatus = []
mesh_010_08e.EdgeMapArrayStatus = []
mesh_010_08e.FaceMapArrayStatus = []
mesh_010_08e.ElementMapArrayStatus = []
mesh_010_08e.ElementBlocks = []
mesh_010_08e.FaceBlocks = []
mesh_010_08e.HasModeShapes = 0
mesh_010_08e.ModeShape = 1
mesh_010_08e.AnimateVibrations = 1
mesh_010_08e.IgnoreFileTime = 0
mesh_010_08e.GenerateFileIdArray = 0
mesh_010_08e.UseLegacyBlockNamesWithElementTypes = 0

# Properties modified on mesh_010_08e
mesh_010_08e.ElementBlocks = ['Unnamed block ID: 0']
mesh_010_08e.FilePrefix = ''
mesh_010_08e.FilePattern = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
mesh_010_08eDisplay = Show(mesh_010_08e, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
mesh_010_08eDisplay.Selection = None
mesh_010_08eDisplay.Representation = 'Surface'
mesh_010_08eDisplay.ColorArrayName = [None, '']
mesh_010_08eDisplay.LookupTable = None
mesh_010_08eDisplay.MapScalars = 1
mesh_010_08eDisplay.MultiComponentsMapping = 0
mesh_010_08eDisplay.InterpolateScalarsBeforeMapping = 1
mesh_010_08eDisplay.Opacity = 1.0
mesh_010_08eDisplay.PointSize = 2.0
mesh_010_08eDisplay.LineWidth = 1.0
mesh_010_08eDisplay.RenderLinesAsTubes = 0
mesh_010_08eDisplay.RenderPointsAsSpheres = 0
mesh_010_08eDisplay.Interpolation = 'Gouraud'
mesh_010_08eDisplay.Specular = 0.0
mesh_010_08eDisplay.SpecularColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.SpecularPower = 100.0
mesh_010_08eDisplay.Luminosity = 0.0
mesh_010_08eDisplay.Ambient = 0.0
mesh_010_08eDisplay.Diffuse = 1.0
mesh_010_08eDisplay.Roughness = 0.3
mesh_010_08eDisplay.Metallic = 0.0
mesh_010_08eDisplay.EdgeTint = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.SelectTCoordArray = 'None'
mesh_010_08eDisplay.SelectNormalArray = 'None'
mesh_010_08eDisplay.SelectTangentArray = 'None'
mesh_010_08eDisplay.Texture = None
mesh_010_08eDisplay.RepeatTextures = 1
mesh_010_08eDisplay.InterpolateTextures = 0
mesh_010_08eDisplay.SeamlessU = 0
mesh_010_08eDisplay.SeamlessV = 0
mesh_010_08eDisplay.UseMipmapTextures = 0
mesh_010_08eDisplay.BaseColorTexture = None
mesh_010_08eDisplay.NormalTexture = None
mesh_010_08eDisplay.NormalScale = 1.0
mesh_010_08eDisplay.MaterialTexture = None
mesh_010_08eDisplay.OcclusionStrength = 1.0
mesh_010_08eDisplay.EmissiveTexture = None
mesh_010_08eDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.FlipTextures = 0
mesh_010_08eDisplay.BackfaceRepresentation = 'Follow Frontface'
mesh_010_08eDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.BackfaceOpacity = 1.0
mesh_010_08eDisplay.Position = [0.0, 0.0, 0.0]
mesh_010_08eDisplay.Scale = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.Orientation = [0.0, 0.0, 0.0]
mesh_010_08eDisplay.Origin = [0.0, 0.0, 0.0]
mesh_010_08eDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
mesh_010_08eDisplay.Pickable = 1
mesh_010_08eDisplay.Triangulate = 0
mesh_010_08eDisplay.UseShaderReplacements = 0
mesh_010_08eDisplay.ShaderReplacements = ''
mesh_010_08eDisplay.NonlinearSubdivisionLevel = 1
mesh_010_08eDisplay.UseDataPartitions = 0
mesh_010_08eDisplay.OSPRayUseScaleArray = 'All Approximate'
mesh_010_08eDisplay.OSPRayScaleArray = 'GlobalNodeId'
mesh_010_08eDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
mesh_010_08eDisplay.OSPRayMaterial = 'None'
mesh_010_08eDisplay.Orient = 0
mesh_010_08eDisplay.OrientationMode = 'Direction'
mesh_010_08eDisplay.SelectOrientationVectors = 'None'
mesh_010_08eDisplay.Scaling = 0
mesh_010_08eDisplay.ScaleMode = 'No Data Scaling Off'
mesh_010_08eDisplay.ScaleFactor = 0.1
mesh_010_08eDisplay.SelectScaleArray = 'GlobalNodeId'
mesh_010_08eDisplay.GlyphType = 'Arrow'
mesh_010_08eDisplay.UseGlyphTable = 0
mesh_010_08eDisplay.GlyphTableIndexArray = 'GlobalNodeId'
mesh_010_08eDisplay.UseCompositeGlyphTable = 0
mesh_010_08eDisplay.UseGlyphCullingAndLOD = 0
mesh_010_08eDisplay.LODValues = []
mesh_010_08eDisplay.ColorByLODIndex = 0
mesh_010_08eDisplay.GaussianRadius = 0.005
mesh_010_08eDisplay.ShaderPreset = 'Sphere'
mesh_010_08eDisplay.CustomTriangleScale = 3
mesh_010_08eDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
mesh_010_08eDisplay.Emissive = 0
mesh_010_08eDisplay.ScaleByArray = 0
mesh_010_08eDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
mesh_010_08eDisplay.ScaleArrayComponent = ''
mesh_010_08eDisplay.UseScaleFunction = 1
mesh_010_08eDisplay.ScaleTransferFunction = 'PiecewiseFunction'
mesh_010_08eDisplay.OpacityByArray = 0
mesh_010_08eDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
mesh_010_08eDisplay.OpacityArrayComponent = ''
mesh_010_08eDisplay.OpacityTransferFunction = 'PiecewiseFunction'
mesh_010_08eDisplay.DataAxesGrid = 'GridAxesRepresentation'
mesh_010_08eDisplay.SelectionCellLabelBold = 0
mesh_010_08eDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
mesh_010_08eDisplay.SelectionCellLabelFontFamily = 'Arial'
mesh_010_08eDisplay.SelectionCellLabelFontFile = ''
mesh_010_08eDisplay.SelectionCellLabelFontSize = 18
mesh_010_08eDisplay.SelectionCellLabelItalic = 0
mesh_010_08eDisplay.SelectionCellLabelJustification = 'Left'
mesh_010_08eDisplay.SelectionCellLabelOpacity = 1.0
mesh_010_08eDisplay.SelectionCellLabelShadow = 0
mesh_010_08eDisplay.SelectionPointLabelBold = 0
mesh_010_08eDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
mesh_010_08eDisplay.SelectionPointLabelFontFamily = 'Arial'
mesh_010_08eDisplay.SelectionPointLabelFontFile = ''
mesh_010_08eDisplay.SelectionPointLabelFontSize = 18
mesh_010_08eDisplay.SelectionPointLabelItalic = 0
mesh_010_08eDisplay.SelectionPointLabelJustification = 'Left'
mesh_010_08eDisplay.SelectionPointLabelOpacity = 1.0
mesh_010_08eDisplay.SelectionPointLabelShadow = 0
mesh_010_08eDisplay.PolarAxes = 'PolarAxesRepresentation'
mesh_010_08eDisplay.ScalarOpacityFunction = None
mesh_010_08eDisplay.ScalarOpacityUnitDistance = 0.0406699310237908
mesh_010_08eDisplay.UseSeparateOpacityArray = 0
mesh_010_08eDisplay.OpacityArrayName = ['POINTS', 'GlobalNodeId']
mesh_010_08eDisplay.OpacityComponent = ''
mesh_010_08eDisplay.ExtractedBlockIndex = 2
mesh_010_08eDisplay.SelectMapper = 'Projected tetra'
mesh_010_08eDisplay.SamplingDimensions = [128, 128, 128]
mesh_010_08eDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
mesh_010_08eDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
mesh_010_08eDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
mesh_010_08eDisplay.GlyphType.TipResolution = 6
mesh_010_08eDisplay.GlyphType.TipRadius = 0.1
mesh_010_08eDisplay.GlyphType.TipLength = 0.35
mesh_010_08eDisplay.GlyphType.ShaftResolution = 6
mesh_010_08eDisplay.GlyphType.ShaftRadius = 0.03
mesh_010_08eDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
mesh_010_08eDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 50413.0, 1.0, 0.5, 0.0]
mesh_010_08eDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
mesh_010_08eDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 50413.0, 1.0, 0.5, 0.0]
mesh_010_08eDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
mesh_010_08eDisplay.DataAxesGrid.XTitle = 'X Axis'
mesh_010_08eDisplay.DataAxesGrid.YTitle = 'Y Axis'
mesh_010_08eDisplay.DataAxesGrid.ZTitle = 'Z Axis'
mesh_010_08eDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.XTitleFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.XTitleBold = 0
mesh_010_08eDisplay.DataAxesGrid.XTitleItalic = 0
mesh_010_08eDisplay.DataAxesGrid.XTitleFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.XTitleShadow = 0
mesh_010_08eDisplay.DataAxesGrid.XTitleOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.YTitleFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.YTitleBold = 0
mesh_010_08eDisplay.DataAxesGrid.YTitleItalic = 0
mesh_010_08eDisplay.DataAxesGrid.YTitleFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.YTitleShadow = 0
mesh_010_08eDisplay.DataAxesGrid.YTitleOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.ZTitleFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.ZTitleBold = 0
mesh_010_08eDisplay.DataAxesGrid.ZTitleItalic = 0
mesh_010_08eDisplay.DataAxesGrid.ZTitleFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.ZTitleShadow = 0
mesh_010_08eDisplay.DataAxesGrid.ZTitleOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.FacesToRender = 63
mesh_010_08eDisplay.DataAxesGrid.CullBackface = 0
mesh_010_08eDisplay.DataAxesGrid.CullFrontface = 1
mesh_010_08eDisplay.DataAxesGrid.ShowGrid = 0
mesh_010_08eDisplay.DataAxesGrid.ShowEdges = 1
mesh_010_08eDisplay.DataAxesGrid.ShowTicks = 1
mesh_010_08eDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
mesh_010_08eDisplay.DataAxesGrid.AxesToLabel = 63
mesh_010_08eDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.XLabelFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.XLabelBold = 0
mesh_010_08eDisplay.DataAxesGrid.XLabelItalic = 0
mesh_010_08eDisplay.DataAxesGrid.XLabelFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.XLabelShadow = 0
mesh_010_08eDisplay.DataAxesGrid.XLabelOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.YLabelFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.YLabelBold = 0
mesh_010_08eDisplay.DataAxesGrid.YLabelItalic = 0
mesh_010_08eDisplay.DataAxesGrid.YLabelFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.YLabelShadow = 0
mesh_010_08eDisplay.DataAxesGrid.YLabelOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
mesh_010_08eDisplay.DataAxesGrid.ZLabelFontFile = ''
mesh_010_08eDisplay.DataAxesGrid.ZLabelBold = 0
mesh_010_08eDisplay.DataAxesGrid.ZLabelItalic = 0
mesh_010_08eDisplay.DataAxesGrid.ZLabelFontSize = 12
mesh_010_08eDisplay.DataAxesGrid.ZLabelShadow = 0
mesh_010_08eDisplay.DataAxesGrid.ZLabelOpacity = 1.0
mesh_010_08eDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
mesh_010_08eDisplay.DataAxesGrid.XAxisPrecision = 2
mesh_010_08eDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
mesh_010_08eDisplay.DataAxesGrid.XAxisLabels = []
mesh_010_08eDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
mesh_010_08eDisplay.DataAxesGrid.YAxisPrecision = 2
mesh_010_08eDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
mesh_010_08eDisplay.DataAxesGrid.YAxisLabels = []
mesh_010_08eDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
mesh_010_08eDisplay.DataAxesGrid.ZAxisPrecision = 2
mesh_010_08eDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
mesh_010_08eDisplay.DataAxesGrid.ZAxisLabels = []
mesh_010_08eDisplay.DataAxesGrid.UseCustomBounds = 0
mesh_010_08eDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
mesh_010_08eDisplay.PolarAxes.Visibility = 0
mesh_010_08eDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
mesh_010_08eDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
mesh_010_08eDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
mesh_010_08eDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
mesh_010_08eDisplay.PolarAxes.EnableCustomRange = 0
mesh_010_08eDisplay.PolarAxes.CustomRange = [0.0, 1.0]
mesh_010_08eDisplay.PolarAxes.PolarAxisVisibility = 1
mesh_010_08eDisplay.PolarAxes.RadialAxesVisibility = 1
mesh_010_08eDisplay.PolarAxes.DrawRadialGridlines = 1
mesh_010_08eDisplay.PolarAxes.PolarArcsVisibility = 1
mesh_010_08eDisplay.PolarAxes.DrawPolarArcsGridlines = 1
mesh_010_08eDisplay.PolarAxes.NumberOfRadialAxes = 0
mesh_010_08eDisplay.PolarAxes.AutoSubdividePolarAxis = 1
mesh_010_08eDisplay.PolarAxes.NumberOfPolarAxis = 0
mesh_010_08eDisplay.PolarAxes.MinimumRadius = 0.0
mesh_010_08eDisplay.PolarAxes.MinimumAngle = 0.0
mesh_010_08eDisplay.PolarAxes.MaximumAngle = 90.0
mesh_010_08eDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
mesh_010_08eDisplay.PolarAxes.Ratio = 1.0
mesh_010_08eDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleVisibility = 1
mesh_010_08eDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
mesh_010_08eDisplay.PolarAxes.PolarLabelVisibility = 1
mesh_010_08eDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
mesh_010_08eDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
mesh_010_08eDisplay.PolarAxes.RadialLabelVisibility = 1
mesh_010_08eDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
mesh_010_08eDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
mesh_010_08eDisplay.PolarAxes.RadialUnitsVisibility = 1
mesh_010_08eDisplay.PolarAxes.ScreenSize = 10.0
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleFontFile = ''
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleBold = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleItalic = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleShadow = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisTitleFontSize = 12
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelFontFile = ''
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelBold = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelItalic = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelShadow = 0
mesh_010_08eDisplay.PolarAxes.PolarAxisLabelFontSize = 12
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextBold = 0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextItalic = 0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextShadow = 0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
mesh_010_08eDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
mesh_010_08eDisplay.PolarAxes.EnableDistanceLOD = 1
mesh_010_08eDisplay.PolarAxes.DistanceLODThreshold = 0.7
mesh_010_08eDisplay.PolarAxes.EnableViewAngleLOD = 1
mesh_010_08eDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
mesh_010_08eDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
mesh_010_08eDisplay.PolarAxes.PolarTicksVisibility = 1
mesh_010_08eDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
mesh_010_08eDisplay.PolarAxes.TickLocation = 'Both'
mesh_010_08eDisplay.PolarAxes.AxisTickVisibility = 1
mesh_010_08eDisplay.PolarAxes.AxisMinorTickVisibility = 0
mesh_010_08eDisplay.PolarAxes.ArcTickVisibility = 1
mesh_010_08eDisplay.PolarAxes.ArcMinorTickVisibility = 0
mesh_010_08eDisplay.PolarAxes.DeltaAngleMajor = 10.0
mesh_010_08eDisplay.PolarAxes.DeltaAngleMinor = 5.0
mesh_010_08eDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
mesh_010_08eDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
mesh_010_08eDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
mesh_010_08eDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
mesh_010_08eDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
mesh_010_08eDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
mesh_010_08eDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
mesh_010_08eDisplay.PolarAxes.ArcMajorTickSize = 0.0
mesh_010_08eDisplay.PolarAxes.ArcTickRatioSize = 0.3
mesh_010_08eDisplay.PolarAxes.ArcMajorTickThickness = 1.0
mesh_010_08eDisplay.PolarAxes.ArcTickRatioThickness = 0.5
mesh_010_08eDisplay.PolarAxes.Use2DMode = 0
mesh_010_08eDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.5, 0.5, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh_010_08eDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
mesh_010_08eDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# change representation type
mesh_010_08eDisplay.SetRepresentationType('Surface With Edges')

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(3200, 2400)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot(outputName, renderView1, ImageResolution=[3200, 2400],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # JPEG options
    Quality=95,
    Progressive=1)
