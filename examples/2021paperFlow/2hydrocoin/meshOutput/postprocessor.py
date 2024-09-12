# trace generated using paraview version 5.9.0-RC4

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'ExodusIIReader'
mesh_040_02e = ExodusIIReader(registrationName='mesh_040_02.e', FileName=['/Users/favinom/projects/real/parrot2/examples/2021paperFlow/complex/meshOutput/mesh_040_02.e'])
mesh_040_02e.GenerateObjectIdCellArray = 1
mesh_040_02e.GenerateGlobalElementIdArray = 1
mesh_040_02e.ElementVariables = []
mesh_040_02e.FaceVariables = []
mesh_040_02e.EdgeVariables = []
mesh_040_02e.SideSetResultArrayStatus = []
mesh_040_02e.NodeSetResultArrayStatus = []
mesh_040_02e.FaceSetResultArrayStatus = []
mesh_040_02e.EdgeSetResultArrayStatus = []
mesh_040_02e.GenerateGlobalNodeIdArray = 1
mesh_040_02e.ElementSetResultArrayStatus = []
mesh_040_02e.PointVariables = []
mesh_040_02e.GlobalVariables = []
mesh_040_02e.ApplyDisplacements = 1
mesh_040_02e.DisplacementMagnitude = 1.0
mesh_040_02e.EdgeBlocks = []
mesh_040_02e.NodeSetArrayStatus = []
mesh_040_02e.SideSetArrayStatus = []
mesh_040_02e.FaceSetArrayStatus = []
mesh_040_02e.EdgeSetArrayStatus = []
mesh_040_02e.ElementSetArrayStatus = []
mesh_040_02e.NodeMapArrayStatus = []
mesh_040_02e.EdgeMapArrayStatus = []
mesh_040_02e.FaceMapArrayStatus = []
mesh_040_02e.ElementMapArrayStatus = []
mesh_040_02e.ElementBlocks = []
mesh_040_02e.FaceBlocks = []
mesh_040_02e.HasModeShapes = 0
mesh_040_02e.ModeShape = 1
mesh_040_02e.AnimateVibrations = 1
mesh_040_02e.IgnoreFileTime = 0
mesh_040_02e.GenerateFileIdArray = 0
mesh_040_02e.UseLegacyBlockNamesWithElementTypes = 0

# Properties modified on mesh_040_02e
mesh_040_02e.ElementBlocks = ['Unnamed block ID: 0']
mesh_040_02e.FilePrefix = ''
mesh_040_02e.FilePattern = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
mesh_040_02eDisplay = Show(mesh_040_02e, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
mesh_040_02eDisplay.Selection = None
mesh_040_02eDisplay.Representation = 'Surface'
mesh_040_02eDisplay.ColorArrayName = [None, '']
mesh_040_02eDisplay.LookupTable = None
mesh_040_02eDisplay.MapScalars = 1
mesh_040_02eDisplay.MultiComponentsMapping = 0
mesh_040_02eDisplay.InterpolateScalarsBeforeMapping = 1
mesh_040_02eDisplay.Opacity = 1.0
mesh_040_02eDisplay.PointSize = 2.0
mesh_040_02eDisplay.LineWidth = 1.0
mesh_040_02eDisplay.RenderLinesAsTubes = 0
mesh_040_02eDisplay.RenderPointsAsSpheres = 0
mesh_040_02eDisplay.Interpolation = 'Gouraud'
mesh_040_02eDisplay.Specular = 0.0
mesh_040_02eDisplay.SpecularColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.SpecularPower = 100.0
mesh_040_02eDisplay.Luminosity = 0.0
mesh_040_02eDisplay.Ambient = 0.0
mesh_040_02eDisplay.Diffuse = 1.0
mesh_040_02eDisplay.Roughness = 0.3
mesh_040_02eDisplay.Metallic = 0.0
mesh_040_02eDisplay.EdgeTint = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.SelectTCoordArray = 'None'
mesh_040_02eDisplay.SelectNormalArray = 'None'
mesh_040_02eDisplay.SelectTangentArray = 'None'
mesh_040_02eDisplay.Texture = None
mesh_040_02eDisplay.RepeatTextures = 1
mesh_040_02eDisplay.InterpolateTextures = 0
mesh_040_02eDisplay.SeamlessU = 0
mesh_040_02eDisplay.SeamlessV = 0
mesh_040_02eDisplay.UseMipmapTextures = 0
mesh_040_02eDisplay.BaseColorTexture = None
mesh_040_02eDisplay.NormalTexture = None
mesh_040_02eDisplay.NormalScale = 1.0
mesh_040_02eDisplay.MaterialTexture = None
mesh_040_02eDisplay.OcclusionStrength = 1.0
mesh_040_02eDisplay.EmissiveTexture = None
mesh_040_02eDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.FlipTextures = 0
mesh_040_02eDisplay.BackfaceRepresentation = 'Follow Frontface'
mesh_040_02eDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.BackfaceOpacity = 1.0
mesh_040_02eDisplay.Position = [0.0, 0.0, 0.0]
mesh_040_02eDisplay.Scale = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.Orientation = [0.0, 0.0, 0.0]
mesh_040_02eDisplay.Origin = [0.0, 0.0, 0.0]
mesh_040_02eDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
mesh_040_02eDisplay.Pickable = 1
mesh_040_02eDisplay.Triangulate = 0
mesh_040_02eDisplay.UseShaderReplacements = 0
mesh_040_02eDisplay.ShaderReplacements = ''
mesh_040_02eDisplay.NonlinearSubdivisionLevel = 1
mesh_040_02eDisplay.UseDataPartitions = 0
mesh_040_02eDisplay.OSPRayUseScaleArray = 'All Approximate'
mesh_040_02eDisplay.OSPRayScaleArray = 'GlobalNodeId'
mesh_040_02eDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
mesh_040_02eDisplay.OSPRayMaterial = 'None'
mesh_040_02eDisplay.Orient = 0
mesh_040_02eDisplay.OrientationMode = 'Direction'
mesh_040_02eDisplay.SelectOrientationVectors = 'None'
mesh_040_02eDisplay.Scaling = 0
mesh_040_02eDisplay.ScaleMode = 'No Data Scaling Off'
mesh_040_02eDisplay.ScaleFactor = 0.1
mesh_040_02eDisplay.SelectScaleArray = 'GlobalNodeId'
mesh_040_02eDisplay.GlyphType = 'Arrow'
mesh_040_02eDisplay.UseGlyphTable = 0
mesh_040_02eDisplay.GlyphTableIndexArray = 'GlobalNodeId'
mesh_040_02eDisplay.UseCompositeGlyphTable = 0
mesh_040_02eDisplay.UseGlyphCullingAndLOD = 0
mesh_040_02eDisplay.LODValues = []
mesh_040_02eDisplay.ColorByLODIndex = 0
mesh_040_02eDisplay.GaussianRadius = 0.005
mesh_040_02eDisplay.ShaderPreset = 'Sphere'
mesh_040_02eDisplay.CustomTriangleScale = 3
mesh_040_02eDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
mesh_040_02eDisplay.Emissive = 0
mesh_040_02eDisplay.ScaleByArray = 0
mesh_040_02eDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
mesh_040_02eDisplay.ScaleArrayComponent = ''
mesh_040_02eDisplay.UseScaleFunction = 1
mesh_040_02eDisplay.ScaleTransferFunction = 'PiecewiseFunction'
mesh_040_02eDisplay.OpacityByArray = 0
mesh_040_02eDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
mesh_040_02eDisplay.OpacityArrayComponent = ''
mesh_040_02eDisplay.OpacityTransferFunction = 'PiecewiseFunction'
mesh_040_02eDisplay.DataAxesGrid = 'GridAxesRepresentation'
mesh_040_02eDisplay.SelectionCellLabelBold = 0
mesh_040_02eDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
mesh_040_02eDisplay.SelectionCellLabelFontFamily = 'Arial'
mesh_040_02eDisplay.SelectionCellLabelFontFile = ''
mesh_040_02eDisplay.SelectionCellLabelFontSize = 18
mesh_040_02eDisplay.SelectionCellLabelItalic = 0
mesh_040_02eDisplay.SelectionCellLabelJustification = 'Left'
mesh_040_02eDisplay.SelectionCellLabelOpacity = 1.0
mesh_040_02eDisplay.SelectionCellLabelShadow = 0
mesh_040_02eDisplay.SelectionPointLabelBold = 0
mesh_040_02eDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
mesh_040_02eDisplay.SelectionPointLabelFontFamily = 'Arial'
mesh_040_02eDisplay.SelectionPointLabelFontFile = ''
mesh_040_02eDisplay.SelectionPointLabelFontSize = 18
mesh_040_02eDisplay.SelectionPointLabelItalic = 0
mesh_040_02eDisplay.SelectionPointLabelJustification = 'Left'
mesh_040_02eDisplay.SelectionPointLabelOpacity = 1.0
mesh_040_02eDisplay.SelectionPointLabelShadow = 0
mesh_040_02eDisplay.PolarAxes = 'PolarAxesRepresentation'
mesh_040_02eDisplay.ScalarOpacityFunction = None
mesh_040_02eDisplay.ScalarOpacityUnitDistance = 0.08906761047582233
mesh_040_02eDisplay.UseSeparateOpacityArray = 0
mesh_040_02eDisplay.OpacityArrayName = ['POINTS', 'GlobalNodeId']
mesh_040_02eDisplay.OpacityComponent = ''
mesh_040_02eDisplay.ExtractedBlockIndex = 2
mesh_040_02eDisplay.SelectMapper = 'Projected tetra'
mesh_040_02eDisplay.SamplingDimensions = [128, 128, 128]
mesh_040_02eDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
mesh_040_02eDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
mesh_040_02eDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
mesh_040_02eDisplay.GlyphType.TipResolution = 6
mesh_040_02eDisplay.GlyphType.TipRadius = 0.1
mesh_040_02eDisplay.GlyphType.TipLength = 0.35
mesh_040_02eDisplay.GlyphType.ShaftResolution = 6
mesh_040_02eDisplay.GlyphType.ShaftRadius = 0.03
mesh_040_02eDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
mesh_040_02eDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 4687.0, 1.0, 0.5, 0.0]
mesh_040_02eDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
mesh_040_02eDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 4687.0, 1.0, 0.5, 0.0]
mesh_040_02eDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
mesh_040_02eDisplay.DataAxesGrid.XTitle = 'X Axis'
mesh_040_02eDisplay.DataAxesGrid.YTitle = 'Y Axis'
mesh_040_02eDisplay.DataAxesGrid.ZTitle = 'Z Axis'
mesh_040_02eDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.XTitleFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.XTitleBold = 0
mesh_040_02eDisplay.DataAxesGrid.XTitleItalic = 0
mesh_040_02eDisplay.DataAxesGrid.XTitleFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.XTitleShadow = 0
mesh_040_02eDisplay.DataAxesGrid.XTitleOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.YTitleFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.YTitleBold = 0
mesh_040_02eDisplay.DataAxesGrid.YTitleItalic = 0
mesh_040_02eDisplay.DataAxesGrid.YTitleFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.YTitleShadow = 0
mesh_040_02eDisplay.DataAxesGrid.YTitleOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.ZTitleFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.ZTitleBold = 0
mesh_040_02eDisplay.DataAxesGrid.ZTitleItalic = 0
mesh_040_02eDisplay.DataAxesGrid.ZTitleFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.ZTitleShadow = 0
mesh_040_02eDisplay.DataAxesGrid.ZTitleOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.FacesToRender = 63
mesh_040_02eDisplay.DataAxesGrid.CullBackface = 0
mesh_040_02eDisplay.DataAxesGrid.CullFrontface = 1
mesh_040_02eDisplay.DataAxesGrid.ShowGrid = 0
mesh_040_02eDisplay.DataAxesGrid.ShowEdges = 1
mesh_040_02eDisplay.DataAxesGrid.ShowTicks = 1
mesh_040_02eDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
mesh_040_02eDisplay.DataAxesGrid.AxesToLabel = 63
mesh_040_02eDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.XLabelFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.XLabelBold = 0
mesh_040_02eDisplay.DataAxesGrid.XLabelItalic = 0
mesh_040_02eDisplay.DataAxesGrid.XLabelFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.XLabelShadow = 0
mesh_040_02eDisplay.DataAxesGrid.XLabelOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.YLabelFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.YLabelBold = 0
mesh_040_02eDisplay.DataAxesGrid.YLabelItalic = 0
mesh_040_02eDisplay.DataAxesGrid.YLabelFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.YLabelShadow = 0
mesh_040_02eDisplay.DataAxesGrid.YLabelOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
mesh_040_02eDisplay.DataAxesGrid.ZLabelFontFile = ''
mesh_040_02eDisplay.DataAxesGrid.ZLabelBold = 0
mesh_040_02eDisplay.DataAxesGrid.ZLabelItalic = 0
mesh_040_02eDisplay.DataAxesGrid.ZLabelFontSize = 12
mesh_040_02eDisplay.DataAxesGrid.ZLabelShadow = 0
mesh_040_02eDisplay.DataAxesGrid.ZLabelOpacity = 1.0
mesh_040_02eDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
mesh_040_02eDisplay.DataAxesGrid.XAxisPrecision = 2
mesh_040_02eDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
mesh_040_02eDisplay.DataAxesGrid.XAxisLabels = []
mesh_040_02eDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
mesh_040_02eDisplay.DataAxesGrid.YAxisPrecision = 2
mesh_040_02eDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
mesh_040_02eDisplay.DataAxesGrid.YAxisLabels = []
mesh_040_02eDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
mesh_040_02eDisplay.DataAxesGrid.ZAxisPrecision = 2
mesh_040_02eDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
mesh_040_02eDisplay.DataAxesGrid.ZAxisLabels = []
mesh_040_02eDisplay.DataAxesGrid.UseCustomBounds = 0
mesh_040_02eDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
mesh_040_02eDisplay.PolarAxes.Visibility = 0
mesh_040_02eDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
mesh_040_02eDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
mesh_040_02eDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
mesh_040_02eDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
mesh_040_02eDisplay.PolarAxes.EnableCustomRange = 0
mesh_040_02eDisplay.PolarAxes.CustomRange = [0.0, 1.0]
mesh_040_02eDisplay.PolarAxes.PolarAxisVisibility = 1
mesh_040_02eDisplay.PolarAxes.RadialAxesVisibility = 1
mesh_040_02eDisplay.PolarAxes.DrawRadialGridlines = 1
mesh_040_02eDisplay.PolarAxes.PolarArcsVisibility = 1
mesh_040_02eDisplay.PolarAxes.DrawPolarArcsGridlines = 1
mesh_040_02eDisplay.PolarAxes.NumberOfRadialAxes = 0
mesh_040_02eDisplay.PolarAxes.AutoSubdividePolarAxis = 1
mesh_040_02eDisplay.PolarAxes.NumberOfPolarAxis = 0
mesh_040_02eDisplay.PolarAxes.MinimumRadius = 0.0
mesh_040_02eDisplay.PolarAxes.MinimumAngle = 0.0
mesh_040_02eDisplay.PolarAxes.MaximumAngle = 90.0
mesh_040_02eDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
mesh_040_02eDisplay.PolarAxes.Ratio = 1.0
mesh_040_02eDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleVisibility = 1
mesh_040_02eDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
mesh_040_02eDisplay.PolarAxes.PolarLabelVisibility = 1
mesh_040_02eDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
mesh_040_02eDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
mesh_040_02eDisplay.PolarAxes.RadialLabelVisibility = 1
mesh_040_02eDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
mesh_040_02eDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
mesh_040_02eDisplay.PolarAxes.RadialUnitsVisibility = 1
mesh_040_02eDisplay.PolarAxes.ScreenSize = 10.0
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleFontFile = ''
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleBold = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleItalic = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleShadow = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisTitleFontSize = 12
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelFontFile = ''
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelBold = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelItalic = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelShadow = 0
mesh_040_02eDisplay.PolarAxes.PolarAxisLabelFontSize = 12
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextBold = 0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextItalic = 0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextShadow = 0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
mesh_040_02eDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
mesh_040_02eDisplay.PolarAxes.EnableDistanceLOD = 1
mesh_040_02eDisplay.PolarAxes.DistanceLODThreshold = 0.7
mesh_040_02eDisplay.PolarAxes.EnableViewAngleLOD = 1
mesh_040_02eDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
mesh_040_02eDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
mesh_040_02eDisplay.PolarAxes.PolarTicksVisibility = 1
mesh_040_02eDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
mesh_040_02eDisplay.PolarAxes.TickLocation = 'Both'
mesh_040_02eDisplay.PolarAxes.AxisTickVisibility = 1
mesh_040_02eDisplay.PolarAxes.AxisMinorTickVisibility = 0
mesh_040_02eDisplay.PolarAxes.ArcTickVisibility = 1
mesh_040_02eDisplay.PolarAxes.ArcMinorTickVisibility = 0
mesh_040_02eDisplay.PolarAxes.DeltaAngleMajor = 10.0
mesh_040_02eDisplay.PolarAxes.DeltaAngleMinor = 5.0
mesh_040_02eDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
mesh_040_02eDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
mesh_040_02eDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
mesh_040_02eDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
mesh_040_02eDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
mesh_040_02eDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
mesh_040_02eDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
mesh_040_02eDisplay.PolarAxes.ArcMajorTickSize = 0.0
mesh_040_02eDisplay.PolarAxes.ArcTickRatioSize = 0.3
mesh_040_02eDisplay.PolarAxes.ArcMajorTickThickness = 1.0
mesh_040_02eDisplay.PolarAxes.ArcTickRatioThickness = 0.5
mesh_040_02eDisplay.PolarAxes.Use2DMode = 0
mesh_040_02eDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.5, 0.5, 10000.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh_040_02eDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
mesh_040_02eDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# change representation type
mesh_040_02eDisplay.SetRepresentationType('Surface With Edges')

# hide color bar/color legend
mesh_040_02eDisplay.SetScalarBarVisibility(renderView1, False)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1982, 1256)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot('/Users/favinom/projects/real/parrot2/examples/2021paperFlow/complex/meshOutput/image_040_02.png', renderView1, ImageResolution=[1600, 1200],
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
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).