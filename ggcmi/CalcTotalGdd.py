import sys, os
sys.path.append(r"D:\UserData\hoek008\GGCMI\PySrc\pcse")

# import pcse;
from cropinforeader import CropInfoProvider;
#from pcse.pcse.models import Wofost71_PP;
from pcse.base_classes import WeatherDataProvider
from pcse.engine import Engine
from pcse.fileinput.hdf5reader import Hdf5WeatherDataProvider
#from AgmerraNetcdfReader import AgmerraWeatherDataProvider;
from numpy import arange, mean, std;  
from datetime import datetime, date;
import sqlite3;
import time
import logging
# import pdb;

def main():
    # Initialise
    datadir = r"D:/UserData/hoek008/GGCMI/phase1/data/"
    db_location = os.path.join(datadir, "ggcmi.db");
    conn = None;
    try:
        conn = sqlite3.connect(db_location);
    except Exception as e:
        print "SQLite: " + str(e);
    if conn == None: exit(1);
    remotedir = r"Z:/projects/ggcmi/data/AgMERRA"
    fn = r"AgMERRA_1980-01-01_2010-12-31f.hf5"

    try:
        # Loop over the crops
        watersupplies = ["rf", "ir"];
        cropdict = CropInfoProvider.getCrops();
        for crop in cropdict:
            # Loop over the 2 water supply options
            for watersupply in watersupplies:
                # Get the crop name and insert a new record into the crop table
                cropname = cropdict[crop];
                crop_no = insert_crop(conn, cropname, watersupply);
                
                # Get a crop info provider'
                t1 = time.time()
                cip = CropInfoProvider(cropname, watersupply, datadir);
                nvlp = cip.getExtent();
                cropdata = cip.getCropData();
                msg = ("Retrieving data for crop %s, %s took %6.1f seconds" % (cropname, watersupply, time.time()-t1))
                print msg;
                
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
                            wdp = Hdf5WeatherDataProvider(fn, lat, lon, remotedir);
                            
                            # Loop over the years
                            t2 = time.time()
                            tsums = [];
                            for year in get_available_years(wdp):
                                try:
                                    # Get timer data for the current year
                                    timerdata = cip.getTimerData(start_day, end_day, year);
                                    sitedata = {};
                                    soildata = {"SMFCF":0.4};
                                
                                    # Run simulation
                                    pheno = Engine(sitedata, timerdata, soildata, cropdata, wdp, config="Wofost71_PhenoOnly.conf")
                                    pheno.run(days=366)
                                    
                                    # Retrieve result
                                    results = pheno.get_output()
                                    if (results != None) and isinstance(results, list):
                                        print results[-1];
                                    sumresults = pheno.get_summary_output();
                                    print sumresults;
                                    if isinstance(sumresults[0], dict) and isinstance(sumresults[0]["TSUM"], float):
                                        tsums.append(sumresults[0]["TSUM"]);
                                except Exception, e:
                                    print str(e);
                                # end try  
                            # end year 
                            
                            # Insert average etc. into database
                            if len(tsums) > 0:
                                insert_tsum(conn, 1, crop_no, lat, lon, mean(tsums), std(tsums, ddof=1), min(tsums), max(tsums), len(tsums));
                            msg = ("Simulating for crop %s, %s took %6.1f seconds" % (cropname, watersupply, time.time()-t2)) 
                            print msg;
                            exit();
                            
                        except Exception, e:
                            print str(e);                                    
                        finally:
                            if (wdp != None): 
                                wdp.close();
                        # end try        
                    # end lon            
                # end lat     
                cip.close();
            # end watersupply
        # end crop        
    finally:
        conn.close();

def get_available_years(wdp):
    result = []
    if isinstance(wdp, WeatherDataProvider):
        tmpList = [wdp.first_date.year, wdp.last_date.year];
        if wdp.first_date != datetime(wdp.first_date.year, 1, 1).date():
            tmpList[0] = wdp.first_date.year + 1;
        if wdp.last_date != datetime(wdp.last_date.year, 12, 31).date():
            tmpList[1] = wdp.last_date.year - 1;
        result = range(tmpList[0], tmpList[1]);
    return result;

def insert_crop(conn, crop_name, mgmt_code):
    # TODO: this is prone to SQL injection - insecure!
    strSql = "INSERT INTO 'crop' VALUES (NULL, '%s', '%s')" % (crop_name, mgmt_code);
    cursor = conn.cursor();
    cursor.execute(strSql)
    conn.commit();
    result = cursor.lastrowid;
    del cursor;
    return result;
    
def insert_tsum(conn, dataset_id, crop_no, lat, lon, avg, stdev, minval, maxval, numobs):
    # TODO: this is prone to SQL injection - insecure!
    strSql = "INSERT INTO 'tsum' VALUES (%i, %02i, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %02i)";
    strSql = strSql % (dataset_id, crop_no, lat, lon, avg, stdev, minval, maxval, numobs);
    cursor = conn.cursor();
    cursor.execute(strSql)
    conn.commit();
    del cursor;
        
if __name__ == '__main__':
    main();