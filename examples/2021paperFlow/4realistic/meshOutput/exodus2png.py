# trace generated using paraview version 5.9.0-RC4

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

inputName= sys.argv[1]
outputName=sys.argv[2]

# trace generated using paraview version 5.9.0-RC4

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'ExodusIIReader'
mesh_00e = ExodusIIReader(registrationName='mesh_00.e', FileName=inputName)
mesh_00e.GenerateObjectIdCellArray = 1
mesh_00e.GenerateGlobalElementIdArray = 1
mesh_00e.ElementVariables = []
mesh_00e.FaceVariables = []
mesh_00e.EdgeVariables = []
mesh_00e.SideSetResultArrayStatus = []
mesh_00e.NodeSetResultArrayStatus = []
mesh_00e.FaceSetResultArrayStatus = []
mesh_00e.EdgeSetResultArrayStatus = []
mesh_00e.GenerateGlobalNodeIdArray = 1
mesh_00e.ElementSetResultArrayStatus = []
mesh_00e.PointVariables = []
mesh_00e.GlobalVariables = []
mesh_00e.ApplyDisplacements = 1
mesh_00e.DisplacementMagnitude = 1.0
mesh_00e.EdgeBlocks = []
mesh_00e.NodeSetArrayStatus = []
mesh_00e.SideSetArrayStatus = []
mesh_00e.FaceSetArrayStatus = []
mesh_00e.EdgeSetArrayStatus = []
mesh_00e.ElementSetArrayStatus = []
mesh_00e.NodeMapArrayStatus = []
mesh_00e.EdgeMapArrayStatus = []
mesh_00e.FaceMapArrayStatus = []
mesh_00e.ElementMapArrayStatus = []
mesh_00e.ElementBlocks = []
mesh_00e.FaceBlocks = []
mesh_00e.HasModeShapes = 0
mesh_00e.ModeShape = 1
mesh_00e.AnimateVibrations = 1
mesh_00e.IgnoreFileTime = 0
mesh_00e.GenerateFileIdArray = 0
mesh_00e.UseLegacyBlockNamesWithElementTypes = 0

# Properties modified on mesh_00e
mesh_00e.ElementBlocks = ['Unnamed block ID: 0']
mesh_00e.FilePrefix = ''
mesh_00e.FilePattern = ''

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
mesh_00eDisplay = Show(mesh_00e, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
mesh_00eDisplay.Selection = None
mesh_00eDisplay.Representation = 'Surface'
mesh_00eDisplay.ColorArrayName = [None, '']
mesh_00eDisplay.LookupTable = None
mesh_00eDisplay.MapScalars = 1
mesh_00eDisplay.MultiComponentsMapping = 0
mesh_00eDisplay.InterpolateScalarsBeforeMapping = 1
mesh_00eDisplay.Opacity = 1.0
mesh_00eDisplay.PointSize = 2.0
mesh_00eDisplay.LineWidth = 1.0
mesh_00eDisplay.RenderLinesAsTubes = 0
mesh_00eDisplay.RenderPointsAsSpheres = 0
mesh_00eDisplay.Interpolation = 'Gouraud'
mesh_00eDisplay.Specular = 0.0
mesh_00eDisplay.SpecularColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.SpecularPower = 100.0
mesh_00eDisplay.Luminosity = 0.0
mesh_00eDisplay.Ambient = 0.0
mesh_00eDisplay.Diffuse = 1.0
mesh_00eDisplay.Roughness = 0.3
mesh_00eDisplay.Metallic = 0.0
mesh_00eDisplay.EdgeTint = [1.0, 1.0, 1.0]
mesh_00eDisplay.SelectTCoordArray = 'None'
mesh_00eDisplay.SelectNormalArray = 'None'
mesh_00eDisplay.SelectTangentArray = 'None'
mesh_00eDisplay.Texture = None
mesh_00eDisplay.RepeatTextures = 1
mesh_00eDisplay.InterpolateTextures = 0
mesh_00eDisplay.SeamlessU = 0
mesh_00eDisplay.SeamlessV = 0
mesh_00eDisplay.UseMipmapTextures = 0
mesh_00eDisplay.BaseColorTexture = None
mesh_00eDisplay.NormalTexture = None
mesh_00eDisplay.NormalScale = 1.0
mesh_00eDisplay.MaterialTexture = None
mesh_00eDisplay.OcclusionStrength = 1.0
mesh_00eDisplay.EmissiveTexture = None
mesh_00eDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
mesh_00eDisplay.FlipTextures = 0
mesh_00eDisplay.BackfaceRepresentation = 'Follow Frontface'
mesh_00eDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.BackfaceOpacity = 1.0
mesh_00eDisplay.Position = [0.0, 0.0, 0.0]
mesh_00eDisplay.Scale = [1.0, 1.0, 1.0]
mesh_00eDisplay.Orientation = [0.0, 0.0, 0.0]
mesh_00eDisplay.Origin = [0.0, 0.0, 0.0]
mesh_00eDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
mesh_00eDisplay.Pickable = 1
mesh_00eDisplay.Triangulate = 0
mesh_00eDisplay.UseShaderReplacements = 0
mesh_00eDisplay.ShaderReplacements = ''
mesh_00eDisplay.NonlinearSubdivisionLevel = 1
mesh_00eDisplay.UseDataPartitions = 0
mesh_00eDisplay.OSPRayUseScaleArray = 'All Approximate'
mesh_00eDisplay.OSPRayScaleArray = 'GlobalNodeId'
mesh_00eDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
mesh_00eDisplay.OSPRayMaterial = 'None'
mesh_00eDisplay.Orient = 0
mesh_00eDisplay.OrientationMode = 'Direction'
mesh_00eDisplay.SelectOrientationVectors = 'None'
mesh_00eDisplay.Scaling = 0
mesh_00eDisplay.ScaleMode = 'No Data Scaling Off'
mesh_00eDisplay.ScaleFactor = 70.0
mesh_00eDisplay.SelectScaleArray = 'GlobalNodeId'
mesh_00eDisplay.GlyphType = 'Arrow'
mesh_00eDisplay.UseGlyphTable = 0
mesh_00eDisplay.GlyphTableIndexArray = 'GlobalNodeId'
mesh_00eDisplay.UseCompositeGlyphTable = 0
mesh_00eDisplay.UseGlyphCullingAndLOD = 0
mesh_00eDisplay.LODValues = []
mesh_00eDisplay.ColorByLODIndex = 0
mesh_00eDisplay.GaussianRadius = 3.5
mesh_00eDisplay.ShaderPreset = 'Sphere'
mesh_00eDisplay.CustomTriangleScale = 3
mesh_00eDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
mesh_00eDisplay.Emissive = 0
mesh_00eDisplay.ScaleByArray = 0
mesh_00eDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
mesh_00eDisplay.ScaleArrayComponent = ''
mesh_00eDisplay.UseScaleFunction = 1
mesh_00eDisplay.ScaleTransferFunction = 'PiecewiseFunction'
mesh_00eDisplay.OpacityByArray = 0
mesh_00eDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
mesh_00eDisplay.OpacityArrayComponent = ''
mesh_00eDisplay.OpacityTransferFunction = 'PiecewiseFunction'
mesh_00eDisplay.DataAxesGrid = 'GridAxesRepresentation'
mesh_00eDisplay.SelectionCellLabelBold = 0
mesh_00eDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
mesh_00eDisplay.SelectionCellLabelFontFamily = 'Arial'
mesh_00eDisplay.SelectionCellLabelFontFile = ''
mesh_00eDisplay.SelectionCellLabelFontSize = 18
mesh_00eDisplay.SelectionCellLabelItalic = 0
mesh_00eDisplay.SelectionCellLabelJustification = 'Left'
mesh_00eDisplay.SelectionCellLabelOpacity = 1.0
mesh_00eDisplay.SelectionCellLabelShadow = 0
mesh_00eDisplay.SelectionPointLabelBold = 0
mesh_00eDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
mesh_00eDisplay.SelectionPointLabelFontFamily = 'Arial'
mesh_00eDisplay.SelectionPointLabelFontFile = ''
mesh_00eDisplay.SelectionPointLabelFontSize = 18
mesh_00eDisplay.SelectionPointLabelItalic = 0
mesh_00eDisplay.SelectionPointLabelJustification = 'Left'
mesh_00eDisplay.SelectionPointLabelOpacity = 1.0
mesh_00eDisplay.SelectionPointLabelShadow = 0
mesh_00eDisplay.PolarAxes = 'PolarAxesRepresentation'
mesh_00eDisplay.ScalarOpacityFunction = None
mesh_00eDisplay.ScalarOpacityUnitDistance = 265.23227233760343
mesh_00eDisplay.UseSeparateOpacityArray = 0
mesh_00eDisplay.OpacityArrayName = ['POINTS', 'GlobalNodeId']
mesh_00eDisplay.OpacityComponent = ''
mesh_00eDisplay.ExtractedBlockIndex = 2
mesh_00eDisplay.SelectMapper = 'Projected tetra'
mesh_00eDisplay.SamplingDimensions = [128, 128, 128]
mesh_00eDisplay.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
mesh_00eDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
mesh_00eDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
mesh_00eDisplay.GlyphType.TipResolution = 6
mesh_00eDisplay.GlyphType.TipRadius = 0.1
mesh_00eDisplay.GlyphType.TipLength = 0.35
mesh_00eDisplay.GlyphType.ShaftResolution = 6
mesh_00eDisplay.GlyphType.ShaftRadius = 0.03
mesh_00eDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
mesh_00eDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 56.0, 1.0, 0.5, 0.0]
mesh_00eDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
mesh_00eDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 56.0, 1.0, 0.5, 0.0]
mesh_00eDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
mesh_00eDisplay.DataAxesGrid.XTitle = 'X Axis'
mesh_00eDisplay.DataAxesGrid.YTitle = 'Y Axis'
mesh_00eDisplay.DataAxesGrid.ZTitle = 'Z Axis'
mesh_00eDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.XTitleFontFile = ''
mesh_00eDisplay.DataAxesGrid.XTitleBold = 0
mesh_00eDisplay.DataAxesGrid.XTitleItalic = 0
mesh_00eDisplay.DataAxesGrid.XTitleFontSize = 12
mesh_00eDisplay.DataAxesGrid.XTitleShadow = 0
mesh_00eDisplay.DataAxesGrid.XTitleOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.YTitleFontFile = ''
mesh_00eDisplay.DataAxesGrid.YTitleBold = 0
mesh_00eDisplay.DataAxesGrid.YTitleItalic = 0
mesh_00eDisplay.DataAxesGrid.YTitleFontSize = 12
mesh_00eDisplay.DataAxesGrid.YTitleShadow = 0
mesh_00eDisplay.DataAxesGrid.YTitleOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.ZTitleFontFile = ''
mesh_00eDisplay.DataAxesGrid.ZTitleBold = 0
mesh_00eDisplay.DataAxesGrid.ZTitleItalic = 0
mesh_00eDisplay.DataAxesGrid.ZTitleFontSize = 12
mesh_00eDisplay.DataAxesGrid.ZTitleShadow = 0
mesh_00eDisplay.DataAxesGrid.ZTitleOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.FacesToRender = 63
mesh_00eDisplay.DataAxesGrid.CullBackface = 0
mesh_00eDisplay.DataAxesGrid.CullFrontface = 1
mesh_00eDisplay.DataAxesGrid.ShowGrid = 0
mesh_00eDisplay.DataAxesGrid.ShowEdges = 1
mesh_00eDisplay.DataAxesGrid.ShowTicks = 1
mesh_00eDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
mesh_00eDisplay.DataAxesGrid.AxesToLabel = 63
mesh_00eDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.XLabelFontFile = ''
mesh_00eDisplay.DataAxesGrid.XLabelBold = 0
mesh_00eDisplay.DataAxesGrid.XLabelItalic = 0
mesh_00eDisplay.DataAxesGrid.XLabelFontSize = 12
mesh_00eDisplay.DataAxesGrid.XLabelShadow = 0
mesh_00eDisplay.DataAxesGrid.XLabelOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.YLabelFontFile = ''
mesh_00eDisplay.DataAxesGrid.YLabelBold = 0
mesh_00eDisplay.DataAxesGrid.YLabelItalic = 0
mesh_00eDisplay.DataAxesGrid.YLabelFontSize = 12
mesh_00eDisplay.DataAxesGrid.YLabelShadow = 0
mesh_00eDisplay.DataAxesGrid.YLabelOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
mesh_00eDisplay.DataAxesGrid.ZLabelFontFile = ''
mesh_00eDisplay.DataAxesGrid.ZLabelBold = 0
mesh_00eDisplay.DataAxesGrid.ZLabelItalic = 0
mesh_00eDisplay.DataAxesGrid.ZLabelFontSize = 12
mesh_00eDisplay.DataAxesGrid.ZLabelShadow = 0
mesh_00eDisplay.DataAxesGrid.ZLabelOpacity = 1.0
mesh_00eDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
mesh_00eDisplay.DataAxesGrid.XAxisPrecision = 2
mesh_00eDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
mesh_00eDisplay.DataAxesGrid.XAxisLabels = []
mesh_00eDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
mesh_00eDisplay.DataAxesGrid.YAxisPrecision = 2
mesh_00eDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
mesh_00eDisplay.DataAxesGrid.YAxisLabels = []
mesh_00eDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
mesh_00eDisplay.DataAxesGrid.ZAxisPrecision = 2
mesh_00eDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
mesh_00eDisplay.DataAxesGrid.ZAxisLabels = []
mesh_00eDisplay.DataAxesGrid.UseCustomBounds = 0
mesh_00eDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
mesh_00eDisplay.PolarAxes.Visibility = 0
mesh_00eDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
mesh_00eDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
mesh_00eDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
mesh_00eDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
mesh_00eDisplay.PolarAxes.EnableCustomRange = 0
mesh_00eDisplay.PolarAxes.CustomRange = [0.0, 1.0]
mesh_00eDisplay.PolarAxes.PolarAxisVisibility = 1
mesh_00eDisplay.PolarAxes.RadialAxesVisibility = 1
mesh_00eDisplay.PolarAxes.DrawRadialGridlines = 1
mesh_00eDisplay.PolarAxes.PolarArcsVisibility = 1
mesh_00eDisplay.PolarAxes.DrawPolarArcsGridlines = 1
mesh_00eDisplay.PolarAxes.NumberOfRadialAxes = 0
mesh_00eDisplay.PolarAxes.AutoSubdividePolarAxis = 1
mesh_00eDisplay.PolarAxes.NumberOfPolarAxis = 0
mesh_00eDisplay.PolarAxes.MinimumRadius = 0.0
mesh_00eDisplay.PolarAxes.MinimumAngle = 0.0
mesh_00eDisplay.PolarAxes.MaximumAngle = 90.0
mesh_00eDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
mesh_00eDisplay.PolarAxes.Ratio = 1.0
mesh_00eDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
mesh_00eDisplay.PolarAxes.PolarAxisTitleVisibility = 1
mesh_00eDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
mesh_00eDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
mesh_00eDisplay.PolarAxes.PolarLabelVisibility = 1
mesh_00eDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
mesh_00eDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
mesh_00eDisplay.PolarAxes.RadialLabelVisibility = 1
mesh_00eDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
mesh_00eDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
mesh_00eDisplay.PolarAxes.RadialUnitsVisibility = 1
mesh_00eDisplay.PolarAxes.ScreenSize = 10.0
mesh_00eDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
mesh_00eDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
mesh_00eDisplay.PolarAxes.PolarAxisTitleFontFile = ''
mesh_00eDisplay.PolarAxes.PolarAxisTitleBold = 0
mesh_00eDisplay.PolarAxes.PolarAxisTitleItalic = 0
mesh_00eDisplay.PolarAxes.PolarAxisTitleShadow = 0
mesh_00eDisplay.PolarAxes.PolarAxisTitleFontSize = 12
mesh_00eDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
mesh_00eDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
mesh_00eDisplay.PolarAxes.PolarAxisLabelFontFile = ''
mesh_00eDisplay.PolarAxes.PolarAxisLabelBold = 0
mesh_00eDisplay.PolarAxes.PolarAxisLabelItalic = 0
mesh_00eDisplay.PolarAxes.PolarAxisLabelShadow = 0
mesh_00eDisplay.PolarAxes.PolarAxisLabelFontSize = 12
mesh_00eDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
mesh_00eDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
mesh_00eDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
mesh_00eDisplay.PolarAxes.LastRadialAxisTextBold = 0
mesh_00eDisplay.PolarAxes.LastRadialAxisTextItalic = 0
mesh_00eDisplay.PolarAxes.LastRadialAxisTextShadow = 0
mesh_00eDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
mesh_00eDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
mesh_00eDisplay.PolarAxes.EnableDistanceLOD = 1
mesh_00eDisplay.PolarAxes.DistanceLODThreshold = 0.7
mesh_00eDisplay.PolarAxes.EnableViewAngleLOD = 1
mesh_00eDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
mesh_00eDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
mesh_00eDisplay.PolarAxes.PolarTicksVisibility = 1
mesh_00eDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
mesh_00eDisplay.PolarAxes.TickLocation = 'Both'
mesh_00eDisplay.PolarAxes.AxisTickVisibility = 1
mesh_00eDisplay.PolarAxes.AxisMinorTickVisibility = 0
mesh_00eDisplay.PolarAxes.ArcTickVisibility = 1
mesh_00eDisplay.PolarAxes.ArcMinorTickVisibility = 0
mesh_00eDisplay.PolarAxes.DeltaAngleMajor = 10.0
mesh_00eDisplay.PolarAxes.DeltaAngleMinor = 5.0
mesh_00eDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
mesh_00eDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
mesh_00eDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
mesh_00eDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
mesh_00eDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
mesh_00eDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
mesh_00eDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
mesh_00eDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
mesh_00eDisplay.PolarAxes.ArcMajorTickSize = 0.0
mesh_00eDisplay.PolarAxes.ArcTickRatioSize = 0.3
mesh_00eDisplay.PolarAxes.ArcMajorTickThickness = 1.0
mesh_00eDisplay.PolarAxes.ArcTickRatioThickness = 0.5
mesh_00eDisplay.PolarAxes.Use2DMode = 0
mesh_00eDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(mesh_00eDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
mesh_00eDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# change representation type
mesh_00eDisplay.SetRepresentationType('Surface With Edges')

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1982, 1256)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [350.0, 300.0, 1781.0792195838865]
renderView1.CameraFocalPoint = [350.0, 300.0, 0.0]
renderView1.CameraParallelScale = 460.9772228646444

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
renderView1.CameraPosition = [350.0, 300.0, 1781.0792195838865]
renderView1.CameraFocalPoint = [350.0, 300.0, 0.0]
renderView1.CameraParallelScale = 460.9772228646444

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).