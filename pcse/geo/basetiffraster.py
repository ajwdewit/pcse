from raster import Raster
import os
from libtiff.libtiff_ctypes import TIFFFieldInfo, TIFFDataType, FIELD_CUSTOM, add_tags

class BaseTiffRaster(Raster):
    # Constants
    DATAFILEXT = "tif";
    HEADEREXT = "tfw"; # WORLD_EXT?
    __BYTESPERCELL = 4;

    # Data attributes - assign some dummy values for the mean time
    name = "dummy.tif";
    folder = os.getcwd();
    nodatavalue = -9999.0;
    byteorder = 'II'; # Little endian
    roty = 0.0;
    rotx = 0.0;
    pageNumber = [0, 1]
    modelPixelScale = [0.0, 0.0, 0.0]
    modelTiepoint = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ModelTransformation = 16 * [0.0]
    GeoKeyDirectory = 32 * [0]
    GeoDoubleParams = 1.0e-317
    GeoAsciiParams = "GCS_WGS_1984|"
    GdalMetadata = '<GDALMetadata>\n  <Item name="RepresentationType" sample="0">ATHEMATIC</Item>\n</GDALMetadata>'
        
    def __init__(self, filepath, *datatype):
        pass
    
    def open(self, mode, ncols=1, nrows=1, xll=0.0, yll=0.0, cellsize=1.0, nodatavalue=-9999.0):
        super(BaseTiffRaster, self).open(mode)
        extra_tags = [
            TIFFFieldInfo(297, 2, 2, TIFFDataType.TIFF_SHORT, FIELD_CUSTOM, True, False, "PageNumber"),          
            TIFFFieldInfo(33550, 3, 3, TIFFDataType.TIFF_DOUBLE, FIELD_CUSTOM, True, False, "ModelPixelScaleTag"),
            TIFFFieldInfo(33922, 6, 6, TIFFDataType.TIFF_DOUBLE, FIELD_CUSTOM, True, False, "ModelTiepointTag"),
            TIFFFieldInfo(34264, 16, 16, TIFFDataType.TIFF_DOUBLE, FIELD_CUSTOM, True, False, "ModelTransformationTag"),
            TIFFFieldInfo(34735, 32, 32, TIFFDataType.TIFF_SHORT, FIELD_CUSTOM, True, False, "GeoKeyDirectoryTag"),
            TIFFFieldInfo(34736, -1, -1, TIFFDataType.TIFF_DOUBLE, FIELD_CUSTOM, True, False, "GeoDoubleParamsTag"),
            TIFFFieldInfo(34737, -1, -1, TIFFDataType.TIFF_ASCII, FIELD_CUSTOM, True, False, "GeoAsciiParamsTag"),
            TIFFFieldInfo(42112, -1, -1, TIFFDataType.TIFF_ASCII, FIELD_CUSTOM, True, False, "GDAL_METADATA"),
            TIFFFieldInfo(42113, -1, -1, TIFFDataType.TIFF_ASCII, FIELD_CUSTOM, True, False, "GDAL_NODATA")
        ]
        add_tags(extra_tags)
        
    def set_extra_tags(self):
        # The following is in view of the georeferencing
        self.datafile.SetField("PageNumber", self.pageNumber)
        self.datafile.SetField("ModelPixelScaleTag", self.modelPixelScale) 
        self.datafile.SetField("ModelTiepointTag", self.modelTiepoint)
        self.datafile.SetField("ModelTransformationTag", self.ModelTransformation)
        self.datafile.SetField("GeoKeyDirectoryTag", self.GeoKeyDirectory)
        #self.datafile.SetField("GeoDoubleParamsTag", self.GeoDoubleParams)
        self.datafile.SetField("GeoAsciiParamsTag", str(self.GeoAsciiParams))
        self.datafile.SetField("GDAL_METADATA", str(self.GdalMetadata))
        
    def get_extra_tags(self):        
        # The following is in view of the georeferencing  
        self.pageNumber = self.datafile.GetField("PageNumber")   
        self.modelPixelScale = self.datafile.GetField("ModelPixelScaleTag")
        self.modelTiepoint = self.datafile.GetField("ModelTiepointTag")
        self.ModelTransformation = self.datafile.GetField("ModelTransformationTag")
        self.GeoKeyDirectory = self.datafile.GetField("GeoKeyDirectoryTag")
        #self.GeoDoubleParams = self.datafile.GetField("GeoDoubleParamsTag")
        self.GeoAsciiParams = self.datafile.GetField("GeoAsciiParamsTag")
        self.GdalMetadata = self.datafile.GetField("GDAL_METADATA")
                
    def copy_extra_tags(self, rg):
        # Check input
        if not isinstance(rg, BaseTiffRaster):
            raise TypeError("Input is not a raster!")
        
        # Ok, go ahead
        if rg.pageNumber != None: self.pageNumber = rg.pageNumber
        self.modelPixelScale = rg.modelPixelScale
        self.modelTiepoint = rg.modelTiepoint
        if rg.ModelTransformation != None: self.ModelTransformation = rg.ModelTransformation
        self.GeoKeyDirectory = rg.GeoKeyDirectory
        #self.GeoDoubleParams = rg.GeoDoubleParams
        self.GeoAsciiParams = rg.GeoAsciiParams
        self.GdalMetadata = rg.GdalMetadata
        
                