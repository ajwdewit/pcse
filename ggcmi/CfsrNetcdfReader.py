from pcse.pcse.netcdf4 import NetcdfWeatherDataProvider;

class CfsrWeatherDataProvider(NetcdfWeatherDataProvider):
    # Dummy lambda function
    dummy = lambda x: x;
    
    def __init__(self, latitude, longitude, fpath=None, force_update=False):
        # Add optional variables; conversion to VAP is not done by means of a lambda function
        self.netcdf_variables.append("hur");
        self.pcse_variables.append(("VAP", "hur", self.dummy, "hPa"));
        
        NetcdfWeatherDataProvider.__init__(self, "CFSR", latitude, longitude, fpath, force_update);
        
