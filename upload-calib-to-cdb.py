#!/usr/bin/env python
""" 
Write TOF calibrations from file
This script must, and can, be run only from a micenet computer
Reads calibrations from a file and uploads them to the CDB
For protection against unintentional uploads, the filename and the upload lines are commented out
To upload: Set the appropriate filenames and uncomment the _CALIFILE and the set_detector lines
"""

from cdb import CalibrationSuperMouse

########## Set the server names and calibration service

# cdb write server
_CDB_W_SERVER = 'http://172.16.246.25:8080'
_CALI_SM = CalibrationSuperMouse(_CDB_W_SERVER)

print "API Version: " + _CALI_SM.get_version()
print "Server name: " + _CALI_SM.get_name()

print _CALI_SM

#####################################################################
########## Valid-from date ranges for the calibration

# December 2011 calibrations, valid from date
#_TIMESTAMP = "2011-11-23 12:00:00.0"

# EMR run calibrations, valid from date
# for data with new tof2 positions, and tof1 position modded per survey
#_TIMESTAMP = "2013-08-01 00:00:00.0"

# new trigger calibration, valid from date
_TIMESTAMP = "2015-03-08 00:00:00.0"

#####################################################################
## _TYPE: is the calibration type and can be one of tw, trigger, t0
## _DEVICE: is the trigger station
## _CALIFILE: is the calibration file to read for a given _TYPE and _DEVICE
#####################################################################
########## timewalk calibration
_TYPE = 'tw'
_DEVICE = 'tof1'
# set the tw-calib file name and uncomment below
#_CALIFILE = open('tofcalibTW_july2014.txt', 'r')

# read the file
_DATA = _CALIFILE.read()

#print _DEVICE, _TYPE, _TIMESTAMP, _DATA

### set the calibration - uncomment below to upload
#print _CALI_SM.set_detector(_DEVICE, _TYPE, _TIMESTAMP, _DATA) 

# close the file
_CALIFILE.close()

#####################################################################
########## trigger calibration
_TYPE = 'trigger'
_DEVICE = 'tof1'
# set the trigger-calib file name and uncomment below
#_CALIFILE = open('tofcalibTrigger_trTOF1_july2014.txt', 'r')

# read the file
_DATA = _CALIFILE.read()

#print _DEVICE, _TYPE, _TIMESTAMP, _DATA

# upload to CDB - uncomment below to upload
#print _CALI_SM.set_detector(_DEVICE, _TYPE, _TIMESTAMP, _DATA) 

# close the file
_CALIFILE.close()

#####################################################################
########## t0 calibration
_TYPE = 't0'
_DEVICE = 'tof1'
# set the t0-calib file name and uncomment below
#_CALIFILE = open('tofcalibT0_trTOF1_july2014.txt', 'r')

# read the file
_DATA = _CALIFILE.read()

#print _DEVICE, _TYPE, _TIMESTAMP, _DATA

### set the calibration - uncomment below to upload
#print _CALI_SM.set_detector(_DEVICE, _TYPE, _TIMESTAMP, _DATA) 

# close the file
_CALIFILE.close()

#####################################################################
