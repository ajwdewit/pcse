import pcse;
from CropInfoReader import CropInfoProvider;
#from pcse.pcse.models import Wofost71_PP;
from pcse.pcse.engine import Engine
from AgmerraNetcdfReader import AgmerraWeatherDataProvider;
from numpy import arange;  

def main():
    # Loop over the crops
    watersupplies = ["rf", "ir"];
    datadir = r"D:\UserData\hoek008\GGCMI\phase1\data"
    cropdict = CropInfoProvider.getCrops();
    for crop in cropdict:
        # Loop over the 2 water supply options
        for watersupply in watersupplies:
            cropname = cropdict[crop];
            cip = CropInfoProvider(cropname, watersupply, datadir);
            nvlp = cip.getExtent();
            cropdata = cip.getCropData();
            
            # Make sure we can loop over the grid cells from UL down to LR corner
            x_range = arange(nvlp.getMinX() + 0.5*nvlp.dx, nvlp.getMaxX() + 0.5*nvlp.dy, nvlp.dx);
            y_range = arange(nvlp.getMinY() + 0.5*nvlp.dy, nvlp.getMaxY() + 0.5*nvlp.dx, nvlp.dy);
            if (nvlp.xcoords_sort != 'ASC'): x_range = reversed(x_range);
            if (nvlp.ycoords_sort != 'DESC'): y_range = reversed(y_range);
            
            # Loop but make sure not to waste time on pixels that are not interesting
            for lat in y_range:
                for lon in x_range: 
                    # Check the landmask
                    if not cip.get_landmask(lon, lat):
                        continue;
                    
                    # Check the data on the growing season; start_day eq.to -99 means that area is
                    # usu. only little above crop_specific base temperature and end_day equal to -99
                    # means that hot periods need to be avoided - 
                    start_day, end_day = cip.getSeasonDates(lon, lat);
                    if (start_day == -99) or (end_day == -99):
                        continue;
    
                    # Initialise
                    wdp = None
                    try: 
                        # Retrieve the relevant weather data
                        wdp = AgmerraWeatherDataProvider(lon, lat, datadir);
                        
                        # Loop over the years
                        for year in wdp.available_years:
                            try:
                                # Get timer data for the current year
                                timerdata = cip.getTimerData(start_day, end_day, year);
                                sitedata = None;
                                soildata = None;
                            
                                # Run simulation
                                #wofsim = Wofost71_PP(sitedata, timerdata, soildata, cropdata, wdp);
                                pheno = Engine(sitedata, timerdata, soildata, cropdata, wdp, config="Wofost71_PhenoOnly.conf")
                                pheno.run(days=366)
                                results = pheno.get_output()
                                sumresults = pheno.get_summary_output()
                                
                                # Store result
                                grid_no = wdp.get_grid_no(lon, lat);
                            except Exception, e:
                                print str(e);   
                    except Exception, e:
                        print str(e);                                    
                    finally:
                        if (wdp != None): 
                            wdp.close();
    
            cip.close();
    
if __name__ == '__main__':
    main();