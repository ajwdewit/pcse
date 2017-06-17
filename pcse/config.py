# -*- coding: utf-8 -*-
# CC-BY-SA, Alterra, Wageningen-UR
# Jappe Franke (jappe.franke@wur.nl), Allard de Wit (allard.dewit@wur.nl), June 2016
'''config file for geobis pcse webservice; houses db connection, webserver multicore processing configuration and crop configuration params'''
import os, sys
import platform
import socket


def check_system_user():
    """Returns True if the script is running under a system
    user without a username or HOME folder."""
    
    if platform.system() == "Windows":
        user = os.getenv("USERNAME")
        if user is None:
            return True
    elif platform.system() == "Linux":
        user = os.getenv("USER")
        if user is None:
            return True
    else:
        msg = "Platform unsupported for getting user home."
        raise NotImplementedError(msg)
    return False


#mysql specific connection params:
myhost = "d0116632" if socket.getfqdn().endswith("wurnet.nl") else "localhost"
myuser = "isidore"
mypwd = "isidore"
mydb = "isidore"

#db connection string (sqlalchemy)
dbc = 'mysql+pymysql://%s:%s@%s:3306/%s'%(myuser, mypwd, myhost, mydb)
 
#change mysql executable and cpu cores for deployment on AWS or SurfSara
if platform.system() == "Linux":
    myexec = "/usr/bin/mysql"
    cpucores= 2
else:
    myexec = r"D:\wamp\bin\mysql\mysql5.6.17\bin\mysql.exe"
    cpucores= 12

#running on webserver or not?
ws = check_system_user()

#content type setting for webserver:
ct = "Access-Control-Allow-Origin: *\nContent-type: application/json;charset=utf-8\n"
#ct= "Content-type: text/html\n<pre>"

#browser debugging on or of?
browser_debugging_enabled = False

#pcse model params
CROP_END_TYPE = "maturity"
CROP_START_TYPE = "sowing"
MAX_DURATION = 250

#offline url params for testing
#qs= "cropid=10&lat=25.73&lon=89.238&sowday=20160315&varid=2"
qs= "cropid=10&lat=12.16&lon=35.44&sowday=20170311&varid=-188"
#qs= "cropid=3&lat=25.73&lon=89.23&sowday=20161001&varid=-188"
#qs= "cropid=2&lat=25.73&lon=89.23&sowday=20161001&varid=25"
#qs= "cropid=9&lat=25.73&lon=89.23&sowday=20161001&varid=6"
#override cache for testing
co= True

#sort qs for consistency
qs = qs.split("&")
qs.sort()
qs = "&".join(qs)

# paths for webserver to find python modules
currentdir= os.path.dirname(os.path.abspath(__file__))
updir = os.path.dirname(os.path.abspath(currentdir))
sys.path.append(updir)

# Target folder for NetCDF files
NetCDFdir = os.path.join(updir, "db_update" ,"NCDF")

# Limit weather alerts up to 3 days ahead
walim = 3

# API token for telegram messaging through the BOT API
telegram_bot_token = "332829090:AAGL8ECQ54xokWk93ir9ptSboRZc3U-KwoI"
telegram_channel = "@isidorechannel"