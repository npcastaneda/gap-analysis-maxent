# Extract climate to genepool area (Tomato)
# Python code for ArcGis10

# Name: ExtractByMask_Ex_02.py
# Description: Extracts the cells of a raster that correspond with the areas
#    defined by a mask.
# Requirements: Spatial Analyst Extension
# Author: ESRI

# Import system modules
import arcgisscripting
gp = arcgisscripting.create(9.3)
gp.toolbox = "SA"

# Set Snap raster
#gp.SnapRaster = "G:\\ncastaneda\\gap-analysis-tomato\\gap_tomato_old\\maxent_modeling\\climate_data\\esri_grid\\bio_1"
gp.SnapRaster = "G:\\ncastaneda\\clim\\bio_30s_esri\\bio_1"

# Set the extent environment using the Extent class.
#gp.Extent ="-182.0000, -61.0000, -28.0000, 54.0000"
gp.Extent ="-110.0000, -56.0000, -29.0000, 14.0000"

# Set environment settings
gp.workspace = "G:\\ncastaneda\\clim\\bio_30s_esri"

# Set local variables
##inRaster = "elevation"
##inMaskData = "G:\\ncastaneda\\gap-analysis-tomato\\gap_tomato\\inputs\\genepool_area.shp"
inPolygon = "-110 -56;-110 14;-29 14;-29 -56"
OutFolder = "G:\\ncastaneda\\gap-analysis-tomato\\gap_tomato\\maxent_modeling\\climate_data\\esri_grid"
OutFolderAscii = "G:\\ncastaneda\\gap-analysis-tomato\\gap_tomato\\maxent_modeling\\climate_data\\esri_ascii"

# Check out Spatial Analyst extension license
gp.CheckOutExtension("Spatial")

# Get a list of ESRI GRIDs from the workspace and print
# Get a list of grids in the workspace.
#
rasters = gp.ListRasters("bio*", "GRID")

for raster in rasters:
    # Process: Extract
    #gp.ExtractByMask_sa(raster, inMaskData, OutFolder + "\\" + raster)
    gp.ExtractByPolygon_sa(raster, inPolygon, OutFolder + "\\" + raster, "INSIDE")
    # Process: Convert to ascii
    InRaster = OutFolder + "\\" + raster
    OutAsciiFile = OutFolderAscii + "\\" + raster + ".asc"
    gp.RasterToASCII_conversion(InRaster, OutAsciiFile)
    
    # Data notification
    print "done " + raster
    