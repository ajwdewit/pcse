from pcse.pcse.netcdf4 import NetcdfWeatherDataProvider;

class PrincetonWeatherDataProvider(NetcdfWeatherDataProvider):
    # Dummy lambda function
    dummy = lambda x: x;
    
    def __init__(self, latitude, longitude, fpath=None, force_update=False):
        # Add optional variables; conversion to VAP is not done by means of a lambda function
        self.netcdf_variables.append("hus");
        self.pcse_variables.append(("VAP", "hus", self.dummy, "hPa"));
        
        NetcdfWeatherDataProvider.__init__(self, "Princeton", latitude, longitude, fpath, force_update);
