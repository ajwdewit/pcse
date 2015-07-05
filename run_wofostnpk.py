from __future__ import print_function
import sys
import os
import csv
import datetime as dt
import pcse
from pcse.engine import Engine
from pcse.base_classes import ParameterProvider

start_day = dt.date(1976,1,1)

def write_CSV(outfile, data):

    rows = [d for d in data if d is not None]
    with open(outfile, 'wb') as fp:
        firstrow = rows[0]
        writer = csv.DictWriter(fp, sorted(firstrow.keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
def main():
    end_day = dt.date(1976,12,31)
    timer = {"GRID_NO":1,
             "CROP_NO":1,
             "CAMPAIGNYEAR":1976,
             "START_DATE":start_day,
             "END_DATE":end_day,
             "CROP_START_DATE":start_day,
             "CROP_START_TYPE":"emergence",
             "CROP_END_DATE":end_day,
             "CROP_END_TYPE":"maturity",
             "MAX_DURATION":300}
    soil = pcse.fileinput.CABOFileReader('ec4.soil')
    crop = pcse.fileinput.CABOFileReader('wwh102.crop')
    fertil = pcse.fileinput.CABOFileReader('manage.data')
    site = {"SMLIM":0.3, "IFUNRN":0, "SSMAX":0, "SSI":0.,
            "WAV":50, "NOTINF":0}
    site.update(fertil)
    weather = pcse.fileinput.CABOWeatherDataProvider('NL1', fpath=r"D:\UserData\WOFOST Control Centre\METEO\CABOWE")
    parvalues = ParameterProvider(site, timer, soil, crop)
    pw = Engine(parvalues,  weather, config="Wofost71_NPK.conf")
    pw.run(days=300)
    r = pw.get_output()
    print("TAGP should be 4608.467 = %8.3f" % r[-1]["TAGP"])
    print("TWSO should be 1014.083 = %8.3f" % r[-1]["TWSO"])
    print("TWST should be 2678.155 = %8.3f" % r[-1]["TWST"])

    write_CSV("NPK_reference_results.csv", pw._saved_output)

if __name__ == '__main__':
        main()