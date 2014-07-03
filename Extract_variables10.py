# Extract climate to genepool area (Tomato)
# Python code for ArcGis10

# Name: ExtractByMask_Ex_02.py
# Description: Extracts the cells of a raster that correspond with the areas
#    defined by a mask.
# Requirements: Spatial Analyst Extension
# Author: ESRI

# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *

# Set the extent environment using the Extent class.
arcpy.env.extent = arcpy.Extent(-182.0000, -61.0000, -28.0000, 54.0000)

# Set environment settings
env.workspace = "C:\\Users\\ncp148\\Documents\\_geodata\\bioclim\\bio_2-5m_esri"

# Set local variables
##inRaster = "elevation"
inMaskData = "C:\\Users\\ncp148\\Documents\\CPP_CWR\\_collaboration\\_fontagro\\inputs\\genepool_area.shp"

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Get a list of ESRI GRIDs from the workspace and print
rasterList = arcpy.ListRasters("*", "GRID")
for raster in rasterList:

    # Execute ExtractByMask
    outExtractByMask = ExtractByMask(raster, inMaskData)    

    # Save the output 
    outExtractByMask.save("C:\\Users\\ncp148\\Documents\\CPP_CWR\\_collaboration\\_fontagro\\inputs\\bio_2_5m_fontagro" + "\\" + raster)

    print "done " + raster
    





