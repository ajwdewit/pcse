import os
from pcse.exceptions import PCSEError
from netCDF4 import Dataset
from pcse.fileinput import CABOFileReader
from pcse.geo.floatingpointraster import FloatingPointRaster
from pcse.geo.netcdf4envelope2d import Netcdf4Envelope2D
from datetime import date, timedelta
import numpy as np

import run_settings

class CropInfoProvider():
    regions = run_settings.regions

    # CABO name - filename
    # TODO add more elements for different continents ???
    _crop_info_sources = run_settings.crop_info_sources 

    # File - path to *.flt file with land mask
    landmask_grid = run_settings.landmask_grid

    # Netcdf dataset with planting and harvest dates
    _ds = None
    _cropdata = None
    _envelope = None

    def __init__(self, crop, watersupply="rf", fpath=None):
        # match the given crop and region with a tuple from the list
        _fn = self._crop_info_sources[0][2]
        for cropname, _, filename, _ in self._crop_info_sources:
            # for the mean time, don't worry about the region
            if cropname == crop:
                _fn = filename
                break

        # read parameters from the filename stored in _crop_info_sources
        _fn = os.path.join(fpath, run_settings.crop_input_folder, run_settings.cabofile_folder, _fn)
        self._cropdata = CABOFileReader(_fn)

        # Prepare to read planting and harvest data from the netcdf file
        _fn = crop + "_" + watersupply + run_settings.growing_season_file_suffix
        _fn = os.path.join(fpath, run_settings.crop_input_folder, run_settings.growing_season_folder, _fn)

        try:
            # Open the file
            self._ds = Dataset(_fn, 'r')
            self._envelope = Netcdf4Envelope2D(self._ds)
        except Exception as e:
            fn = os.path.basename(self._ds.filepath())
            raise PCSEError("An error occurred while opening file " + fn + " (" + str(e) + ")")

    @staticmethod
    def getCropGroup(aCropName):
        result = 0
        for cropname, _, _, group_no in CropInfoProvider._crop_info_sources:
            if cropname == aCropName:
                result = group_no
                break
        return result

    @staticmethod
    def getCrops():
        result = {}
        i = 1
        for cropname, _, _, _ in CropInfoProvider.self._ds_crop_info_sources:
            result[i] = cropname
            i = i + 1
        return result

    def getCropData(self):
        return self._cropdata

    def getSeasonDates(self, longitude, latitude):
        eps = 0.001

        # Check that the netCDF file is loaded
        if self._ds is None:
            raise PCSEError("file is no more open.")

        # Check that the given latitude and longitude is within extent of the netCDF4 file
        if not self._envelope.isWithinExtent(longitude, latitude):
            raise PCSEError("Given lat-lon coordinates are beyond the borders of the file.")
        k, i = self._envelope.getColAndRowIndex(longitude, latitude)

        # Check that the found indices are really linked to the given lat-lon
        msg = " coordinate for this index not as expected: "
        assert abs(self._ds.variables["lat"][i] - latitude) < 0.5*self._envelope.dy + eps, "Y" + msg + str(i)
        assert abs(self._ds.variables["lon"][k] - longitude) < 0.5*self._envelope.dx + eps, "X" + msg + str(k)

        # Take into account that the variables are made out of masked arrays
        start_doy = -99
        end_doy = -99
        try:
            arr_elem = self._ds.variables['planting day'][i, k]
            start_doy = int(arr_elem)
            if start_doy < 0: start_doy = -99
            arr_elem = self._ds.variables['harvest day'][i, k]
            end_doy = int(arr_elem)
            if end_doy < 0: end_doy = -99
        except np.ma.core.MaskError:
            pass
        return start_doy, end_doy


    def getTimerData(self, start_doy, end_doy, year):

        try:
            # Prepare the timer data
            result = {}
            if end_doy > start_doy:
                # Assume that start and end time are within the same year
                crop_start_date = date(year,1,1) + timedelta(days=start_doy-1)
                crop_end_date = date(year,1,1) + timedelta(days=end_doy-1)
            else:
                # Not within 1 year
                crop_start_date = date(year,1,1) + timedelta(days=start_doy-1)
                crop_end_date = date(year+1,1,1) + timedelta(days=end_doy-1)
            result['START_DATE'] = crop_start_date
            result['END_DATE'] = crop_end_date
            result['CROP_START_DATE'] = crop_start_date
            result['CROP_END_DATE'] = crop_end_date
            result['CAMPAIGNYEAR'] = year
            result['CROP_START_TYPE'] = 'sowing'
            result['CROP_END_TYPE'] = 'harvest'
            result['MAX_DURATION'] = 365
            return result
        except Exception as e:
            raise PCSEError("An error occurred while preparing the timer data: " + str(e))

    def close(self):
        if self._ds != None:
            self._ds.close()
            self._ds = None

    def getExtent(self):
        return self._envelope

    def _getFullPath(self, fname):
        path = os.path.dirname(self._ds.filepath())
        result = os.path.join(path, fname)
        result = os.path.normpath(result)
        return result

    def _get_value_from_grid(self, longitude, latitude, fpath):
        # Open the file. Get right row and column. Elevations are linked to the cell centres
        r = FloatingPointRaster(fpath, "i")
        if not r.open('r'): raise Exception("Unable to open input file " + r.name)
        k, i = r.getColAndRowIndex(longitude, latitude)
        if (i == r.nrows): i = i - 1

        # Now get hold of the right row, read the wanted value and close
        for _ in range(0, i+1):
            rawline = r.next()
        r.close()
        if (k == r.ncols): k = k - 1
        line = np.frombuffer(rawline, 'f')
        return line[int(k)]

    def get_landmask(self, longitude, latitude):
        fpath = self._getFullPath(self.landmask_grid)
        value = self._get_value_from_grid(longitude, latitude, fpath)
        return (value == 1)


