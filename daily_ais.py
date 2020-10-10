''' daily_ais.py--To be run by cron daily on a radial codar station.

Controlled by the specified .ini configuration file.

The loop files found in loopDir (defined in configuration file) are read in only once;
on subsequent days they are ignored.

The optional -f flag forces all the loop files found to be processed, whether or not they have been before.

Syntax: python3 daily_ais.py [-f] cfgFile

e.g.,   python3 daily_ais.py /Users/codar/Documents/bin/codar_configs/daily_ais_vgpt.ini

REVISION HISTORY

2018-09-12, kpb--Backwards compatibility for Python 2 configparser, matplotlib.

2018-09-14, kpb--Changed format of .png file names for easier use with SeaSonde Archivalist.

2018-09-27, kpb--Changed stationName in all output filenames to be upper case. Again, for easier use with SeaSonde Archivalist.

2018-10-02, kpb--Bugfix. Was reading in partial day of loop data and flagging that day as done. Result was much reduced data in pattern.

2018-10-02, kpb--Added -f flag...

2018-10-05, kpb--Fixed permission problem when opening new yearly file.

2018-10-05, kpb--Changed to skip a corrupted loop file rather than bailing out of program..

2018-10-05, kpb--Modified plot_AIS_pattern() to cope with different lengths of loop1/loop2 patterns.

2018-10-05, kpb--Use try/except instead of testing for shelf file existence.

2018-10-07, kpb--Mysterious bug: after manual "-f" run to clear backlog of loop files, subsequent cron-driven calls to daily_ais.py overwrite the existing pattern-history shelf file. Added numerous diagnostic print() statements to try to solve this.

2018-10-10, kpb--Added test for existence of shelf file before a new, blank one is created. Hopefully will prevent mysterious over-writing of existing shelf file (see above).

2018-10-10, kpb--Refactored. Put calls to plotting routines into generate_plots() function.

2020-10-09, kpb. Added test for corrupt patternHistory shelve file.

'''

class DailyAisError(Exception):
   """Basic exception for daily_ais module."""

import sys
from Codar_loop import Codar_loop
from pattern import Pattern
from pattern import PatternHistory
from pattern import read_pattern_file
import matplotlib.pyplot as plt
from matplotlib import gridspec
from math import radians
from math import degrees
from math import cos
from math import sin
from math import asin
from math import atan2
import numpy as np
import shelve,pickle
from dbm import error as ShelfOpenError
import datetime as dt
import glob
import os
import argparse

try:
  import configparser
except ImportError:
  # Importing configparser failed, presumably because running python2, which
  # has "ConfigParser".
  print("Failed configparser import. Import ConfigParser instead.")
  import ConfigParser as configparser
except Exception as e:
  print("Unexpected error. Failed to import configparser/ConfigParser")
  raise e

scriptName = os.path.basename(__file__)

LOOP1_COLOUR = 'red'
LOOP2_COLOUR = 'blue'

# Useful for debugging, but may cause program to terminate for 
# problems that aren't critical. Comment out for normal use:
#import warnings
#warnings.filterwarnings("error")

################################################################
def Rnb2LL(lat1,lon1,targetRange,targetBearing):
  # Converts a range [km] and a bearing [degrees True] from a
  # latitude/longitude location to a second lat/lon.
  # Algorithm from https://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing
  R = 6378.1 #Radius of the Earth
  lat1 = radians(lat1)
  lon1 = radians(lon1)
  targetBearing = radians(targetBearing)
  lat2 = asin( sin(lat1)*cos(targetRange/R) + cos(lat1)*sin(targetRange/R)*cos(targetBearing))
  lon2 = lon1 + atan2(sin(targetBearing)*sin(targetRange/R)*cos(lat1), cos(targetRange/R)-sin(lat1)*sin(lat2))
  lat2 = degrees(lat2)
  lon2 = degrees(lon2)
  return lat2, lon2


################################################################
def plot_coverage(ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser):
    DOMAP = cfgParser.getboolean('map','DOMAP')

    stationMarkerSize = 40
    #plt.axes(ax)
    plt.sca(ax)

    if DOMAP:
        m.drawcoastlines()
        m.fillcontinents(color='#cc9966',lake_color='#99ffff',zorder=0)

        # Plot location of radar.
        m.scatter(station_x, station_y,s=stationMarkerSize,marker='o',color='b',ax=ax)
        # Plot non-filtered and filtered data.
        m.scatter(ais_x, ais_y,s=1,marker='.',color='gray',ax=ax)
        m.scatter(ais_x[thisDoKeep], ais_y[thisDoKeep],s=1,marker='.',color='red',ax=ax)
    else:
        # Plot non-filtered and filtered data.
        ax.scatter(x=ais_x, y=ais_y,s=1,marker='.',color='gray')
        ax.scatter(x=ais_x[thisDoKeep], y=ais_y[thisDoKeep],s=1,marker='.',color='red')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)

################################################################
def plot_filter_results(loopObj,date,cfgParser):

   filters = get_filters(cfgParser)

   stationLat = cfgParser.getfloat('station','stationLat')
   stationLon = cfgParser.getfloat('station','stationLon')
   stationName = cfgParser.get('station','name')
   mapCentreLat = cfgParser.getfloat('map','mapCentreLat')
   mapCentreLon = cfgParser.getfloat('map','mapCentreLon')

   # Set up variables for map/polar plot.
   targetBearing = np.array(loopObj.get_column('TRGB'))
   targetDistance_m = np.array(loopObj.get_column('TRGD'))
   targetDistance_km = np.array([x/1000. for x in targetDistance_m])

   DOMAP = cfgParser.getboolean('map','DOMAP')

   if DOMAP:
       axProjection = 'rectilinear'
       forceBaseMapSetup = cfgParser.getboolean('map','forceBaseMapSetup')
       mapCacheFile = cfgParser.get('map','mapCacheFile')
       if forceBaseMapSetup or (not os.path.exists(mapCacheFile)):
           # Setting up the basemap must be done if
           #   (a) the config file requires it; or
           #   (b) the basemap cache file does not exist, or
           #   (c) both
           mapWidthMetres = cfgParser.getfloat('map','mapWidthMetres')
           mapHeightMetres = cfgParser.getfloat('map','mapHeightMetres')
           #m = Basemap(projection='tmerc', lat_0=mapCentreLat, lon_0=mapCentreLon, resolution='h', width=75000.0,height=75000.0)
           m = Basemap(projection='tmerc', lat_0=mapCentreLat, lon_0=mapCentreLon, resolution='h', width=mapWidthMetres,height=mapHeightMetres)

           try:
               print("Saving basemap to cache file",mapCacheFile)
               fid = open(mapCacheFile,'wb')
               pickle.dump(m,fid)
               fid.close()
           except:
               #print("Failed to save basemap cache file.")
               raise DailyAisError("Failed to save basemap cache file.")

       else:
           try:
               print("Loading basemap cache from ",mapCacheFile)
               fid = open(mapCacheFile,'rb')
               m = pickle.load(fid)
               fid.close()
               print("loading cache completed")

           except:
               #print("Failed to load pickled map cache file")
               raise DailyAisError("Failed to load pickled map cache file")

       station_x, station_y = m(stationLon,stationLat)

       # Convert AIS target range and bearings to lat, lon.
       aisLats = []
       aisLons = []
       for iTarg in range(len(targetDistance_km)):
           lat2,lon2 = Rnb2LL(stationLat,stationLon,targetDistance_km[iTarg],targetBearing[iTarg])
           aisLats.append(lat2)
           aisLons.append(lon2)
       aisLats = np.array(aisLats)
       aisLons = np.array(aisLons)

       # Convert AIS target locations to map coordinates.
       ais_x, ais_y = m(aisLons,aisLats)

   else:
       axProjection = 'polar'
       m = []
       station_x = []
       station_y = []
       ais_x = np.array([radians(x) for x in targetBearing])
       ais_y = targetDistance_km

   filterTitleNudge = .98
   filterTitleFontSize = 8
   gridX = []
   gridX.extend(filters['bearingStart'])
   gridX.extend(filters['bearingEnd'])
   gridX = np.array(gridX)

   # Set up figure and axes.

   # ...Figure height is too big for a normal monitor, so figure is generated
   # and saved without creating a window for it (in Spyder IDE, can have it
   # appear as an "inline" figure). If plt.show() is used, python will
   # distort the plot to fit it into your computer's screen, which may be
   # undesirable.
   FIG_WIDTH_INCHES = 8
   FIG_HEIGHT_INCHES = 24
   fig = plt.figure(figsize=(FIG_WIDTH_INCHES,FIG_HEIGHT_INCHES))
   gs = gridspec.GridSpec(7, 4)
   master_bearing_ax = fig.add_subplot(gs[0,0:2])
   master_map_ax = fig.add_subplot(gs[0,2:4],projection=axProjection)
   snr_bearing_ax = fig.add_subplot(gs[1,0:2])
   snr_map_ax = fig.add_subplot(gs[1,2:4],projection=axProjection)
   rng_bearing_ax = fig.add_subplot(gs[2,0:2])
   rng_map_ax = fig.add_subplot(gs[2,2:4],projection=axProjection)
   spd_bearing_ax = fig.add_subplot(gs[3,0:2])
   spd_map_ax = fig.add_subplot(gs[3,2:4],projection=axProjection)
   nearBragg_bearing_ax = fig.add_subplot(gs[4,0:2])
   nearBragg_map_ax = fig.add_subplot(gs[5,2:4],projection=axProjection)
   doppWidth_bearing_ax = fig.add_subplot(gs[5,0:2])
   doppWidth_map_ax = fig.add_subplot(gs[4,2:4],projection=axProjection)
   brg_bearing_ax = fig.add_subplot(gs[6,0:2])
   brg_map_ax = fig.add_subplot(gs[6,2:4],projection=axProjection)
   gs.update(wspace=0.1, hspace=0.2)

   # Plot the data.
   titleStr = stationName.upper() + " single-day filtered AIS data, " + str(date)
   fig.suptitle(titleStr, fontsize=14)

   thisDoKeep = loopObj.doKeep['master']
   thisY = np.array(loopObj.get_column('TRGB'))
   #plt.axes(master_bearing_ax)
   plt.sca(master_bearing_ax)

   # Can't plot histogram of unfiltered bearings on same axes as histogram of
   # filtered bearings because the percentage kept is so small that the
   # filtered bearings aren't visible.
   # bearingBins = np.arange(0,360,1)
   #n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)
   #n, bins, patches = plt.hist(thisY, 50, normed=0, facecolor='gray', alpha=0.75)
   #n, bins, patches = plt.hist(thisY[thisDoKeep], np.arange(0,360,5), normed=0, facecolor='red', alpha=0.75)
   n, bins, patches = plt.hist(thisY[thisDoKeep], np.arange(0,360,5), facecolor='red', alpha=0.75)

   master_bearing_ax.set_ylabel('#targets', fontsize=12)
   master_bearing_ax.set_xlim([0, 360])
   master_bearing_ax.set_xticks([0,90,180,270,360],minor=False)
   master_bearing_ax.set_xticks(gridX,minor=True)
   master_bearing_ax.xaxis.grid(True, which='minor')
   master_bearing_ax.set_xticklabels(())
   plot_coverage(master_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)
   #titleStr = stationName.upper() + " single-day filtered AIS data, " + str(date)
   numPts = len(thisY)
   numKept = np.size(np.nonzero(thisDoKeep))
   percentKept = 100 * numKept/numPts
   titleStr = str(numKept) + " points of " + str(len(thisY)) + " (" + "%0.1f" % percentKept + " percent) kept"
   master_bearing_ax.set_title(titleStr,fontsize=filterTitleFontSize,y=filterTitleNudge)

   # The SNR keep index is derived from all three antennas, but don't want to
   # multiply axes unnecessarily. Plot only monopole (but filtered by all 3).
   thisY = np.array(loopObj.get_column('A3SN'))
   thisDoKeep = np.logical_and(loopObj.doKeep['minSNR'],loopObj.doKeep['maxSNR'])
   snr_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   snr_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   snr_bearing_ax.set_xlim([0, 360])
   snr_bearing_ax.set_xticks([0,90,180,270,360])
   snr_bearing_ax.set_xticklabels(())
   snr_bearing_ax.set_ylabel('SNR', fontsize=12)
   snr_bearing_ax.set_title("Antenna 3 SNR, filtered by SNR limits (of all 3 antennas)", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(snr_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)

   thisY = targetDistance_km
   thisDoKeep = np.logical_and(loopObj.doKeep['minRange'],loopObj.doKeep['maxRange'])
   rng_bearing_ax.set_ylabel('RNG [km]', fontsize=12)
   rng_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   rng_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   rng_bearing_ax.set_xlim([0, 360])
   rng_bearing_ax.set_xticks([0,90,180,270,360])
   rng_bearing_ax.set_xticklabels(())
   rng_bearing_ax.set_title("Target distance, filtered by distance limits", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(rng_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)

   thisY = np.absolute(np.array(loopObj.get_column('TRGV')))
   thisDoKeep = np.logical_and(loopObj.doKeep['minVesselSpeed'],loopObj.doKeep['maxVesselSpeed'])
   spd_bearing_ax.set_ylabel('Spd [m/s]', fontsize=12)
   spd_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   spd_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   spd_bearing_ax.set_xlim([0, 360])
   spd_bearing_ax.set_xticks([0,90,180,270,360])
   spd_bearing_ax.set_xticklabels(())
   spd_bearing_ax.set_title("Target speed, filtered by speed limits", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(spd_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)

   thisDoKeep = loopObj.doKeep['nearBragg']
   nearBragg_bearing_ax.set_ylabel('Spd [m/s]', fontsize=12)
   nearBragg_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   nearBragg_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   nearBragg_bearing_ax.set_xlim([0, 360])
   nearBragg_bearing_ax.set_xticks([0,90,180,270,360])
   nearBragg_bearing_ax.set_xticklabels(())
   nearBragg_bearing_ax.set_title("Target speed, filtered by near-Bragg speeds", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(nearBragg_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)

   thisY = np.array(loopObj.get_column('PKD2')) - np.array(loopObj.get_column('PKD1'))
   thisDoKeep = loopObj.doKeep['maxDoppWidth']
   doppWidth_bearing_ax.set_ylabel('doppWidth', fontsize=12)
   doppWidth_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   doppWidth_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   doppWidth_bearing_ax.set_xlim([0, 360])
   doppWidth_bearing_ax.set_xticks([0,90,180,270,360])
   doppWidth_bearing_ax.set_xticklabels(())
   doppWidth_bearing_ax.set_title("Doppler width, filtered by maximum Doppler width", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(doppWidth_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)

   thisY = np.array(loopObj.get_column('TRGB'))
   thisDoKeep = loopObj.doKeep['targetBearing']
   brg_bearing_ax.set_ylabel('brg', fontsize=12)
   brg_bearing_ax.scatter(x=targetBearing, y=thisY, marker='.',s=1,color='gray')
   brg_bearing_ax.scatter(x=targetBearing[thisDoKeep], y=thisY[thisDoKeep], marker='.',s=1,color='red')
   brg_bearing_ax.set_xlim([0, 360])
   brg_bearing_ax.set_xticks(gridX,minor=True)
   brg_bearing_ax.xaxis.grid(True, which='minor')
   brg_bearing_ax.set_xticks([0,90,180,270,360])
   brg_bearing_ax.set_title("Target bearing, filtered by bearing limits", fontsize=filterTitleFontSize,y=filterTitleNudge)
   plot_coverage(brg_map_ax,m,station_x, station_y,ais_x, ais_y,thisDoKeep,cfgParser)
   brg_bearing_ax.set_xlabel('Target Bearing [degrees True]', fontsize=12)

   #fig.savefig('test.png',dpi=200)
   print("plot_filter_results complete")
   return fig

################################################################
def plot_AIS_history(loop1PatternHistory,loop2PatternHistory,loop1MeasPattern,loop2MeasPattern,lastDate,timeAveragePlotLength_days,cfgParser):
    #print("History dates for plotting: ",loop1PatternHistory.patterns.keys())
    loop1PatternHistorySubset = loop1PatternHistory.subset(timeAveragePlotLength_days,lastDate)
    loop2PatternHistorySubset = loop2PatternHistory.subset(timeAveragePlotLength_days,lastDate)
    numDaysAveraged = len(loop1PatternHistorySubset.patterns)

    binSize_degrees = cfgParser.getfloat('processing','binSize_degrees')
    smoothedBinSize_degrees = cfgParser.getfloat('processing','smoothedBinSize_degrees')

    # For loop1, create a single pattern consisting of all the elements in the history.
    bearings, amplitudes, weights = loop1PatternHistorySubset.get_arrays()
    loop1Pattern = Pattern(bearings,amplitudes,binSize_degrees,weights)

    # Create a smoothed pattern for loop 1.
    smoothedLoop1Pattern = Pattern(bearings,amplitudes,smoothedBinSize_degrees,weights)

    # Do the same for loop 2.
    bearings, amplitudes, weights = loop2PatternHistorySubset.get_arrays()
    loop2Pattern = Pattern(bearings,amplitudes,binSize_degrees,weights)
    smoothedLoop2Pattern = Pattern(bearings,amplitudes,smoothedBinSize_degrees,weights)

    # Pass the extracted data to plot_AIS_pattern()
    #fig = plot_AIS_pattern(A13pattern,A23pattern,smoothedA13pattern,smoothedA23pattern, loop1MeasPattern,loop2MeasPattern,thisDay,numDaysAveraged,cfgParser)
    fig = plot_AIS_pattern(loop1Pattern,loop2Pattern,smoothedLoop1Pattern,smoothedLoop2Pattern, loop1MeasPattern,loop2MeasPattern,lastDate,numDaysAveraged,cfgParser)
    return fig

################################################################
def plot_AIS_pattern(A13pattern,A23pattern,smoothedA13pattern,smoothedA23pattern, loop1MeasPattern,loop2MeasPattern,date,numDaysAveraged,cfgParser):
   # 2018-10-05, kpb--Separated out targetBearing1/targetBearing2, smoothedBearing1/smoothedBearing2, etc.
   # Necessitated by new version of average_by_bearing_bins() that removes non-finite data points,
   # sometimes resulting in patterns 1 and 2 having different numbers of points.
   print("starting plot_AIS_pattern")
   stationName = cfgParser.get('station','name')

   targetBearing1, amp1, weights = A13pattern.get_arrays()
   #print("A13 len(targetBearing), len(amp1)",len(targetBearing), len(amp1))
   targetBearing2, amp2, weights = A23pattern.get_arrays()
   #print("A23 len(targetBearing), len(amp2)",len(targetBearing), len(amp2))

   smoothedBearing1, smoothedAmp1, weights = smoothedA13pattern.get_arrays()
   smoothedBearing2, smoothedAmp2, weights = smoothedA23pattern.get_arrays()
   measBearing1, measAmp1, weights = loop1MeasPattern.get_arrays()
   measBearing2, measAmp2, weights = loop2MeasPattern.get_arrays()

   # Where zero (True North) crossings occur, plots are messy unless sorted by increasing bearing.
   sortInd1 = smoothedBearing1.argsort()
   sortInd2 = smoothedBearing2.argsort()

   fig = plt.figure(figsize=(18,6.5))
   gs = gridspec.GridSpec(2, 4)

   ax1 = fig.add_subplot(gs[0,0])

   x = targetBearing1
   y = np.real(amp1)

   ax1.scatter(x=targetBearing1, y=np.real(amp1), marker='.',s=1,color='gray')
   ax1.plot(smoothedBearing1[sortInd1], np.real(smoothedAmp1[sortInd1]), marker='s',markerfacecolor='none',color=LOOP1_COLOUR)
   ax1.set_ylabel('Real(A13)', fontsize=12)
   ax1.set_xlim([0, 360])
   ax1.set_xticks([0,90,180,270,360])
   titleStr = "Amplitude components, averaged to " + str(smoothedA13pattern.binSizeDeg) + "-degree bins"

   try:
     ax1.set_title(titleStr, fontsize=10,loc='left')
   except AttributeError:
     # Python 2 failing with attribute "loc".
     ax1.set_title(titleStr, fontsize=10,horizontalalignment='left')

   ax2 = fig.add_subplot(gs[0,1])
   ax2.scatter(x=targetBearing1, y=np.imag(amp1), marker='.',s=1,color='gray')
   #ax2.scatter(x=bearingBins, y=np.imag(A13_averaged_complex), marker='o',s=2,color=LOOP1_COLOUR)
   ax2.plot(smoothedBearing1[sortInd1], np.imag(smoothedAmp1[sortInd1]), marker='s',markerfacecolor='none',color=LOOP1_COLOUR)
   ax2.set_ylabel('Imag(A13)', fontsize=12)
   ax2.set_xlim([0, 360])
   ax2.set_xticks([0,90,180,270,360])

   ax3 = fig.add_subplot(gs[1,0])
   ax3.scatter(x=targetBearing2, y=np.real(amp2), marker='.',s=1,color='gray')
   #ax3.scatter(x=bearingBins, y=np.real(A23_averaged_complex), marker='o',s=2,color=LOOP2_COLOUR)
   ax3.plot(smoothedBearing2[sortInd2], np.real(smoothedAmp2[sortInd2]), marker='s',markerfacecolor='none',color=LOOP2_COLOUR)
   ax3.set_ylabel('Real(A23)', fontsize=12)
   ax3.set_xlim([0, 360])
   ax3.set_xticks([0,90,180,270,360])

   ax4 = fig.add_subplot(gs[1,1])
   ax4.scatter(x=targetBearing2, y=np.imag(amp2), marker='.',s=1,color='gray')
   #ax4.scatter(x=bearingBins, y=np.imag(A23_averaged_complex), marker='o',s=2,color=LOOP2_COLOUR)
   ax4.plot(smoothedBearing2[sortInd2], np.imag(smoothedAmp2[sortInd2]), marker='s',markerfacecolor='none',color=LOOP2_COLOUR)
   ax4.set_ylabel('Imag(A23)', fontsize=12)
   ax4.set_xlim([0, 360])
   ax4.set_xticks([0,90,180,270,360])

   # Plot the antenna patterns on polar axes.
   ax5 = fig.add_subplot(gs[0:,2:4], polar=True)
   ax5.set_theta_zero_location('N')
   ax5.set_theta_direction(-1)
   ax5.set_yticklabels([])

   try:
     ax5.set_title("Antenna patterns: in use (lines) and AIS-measured (circles)", fontsize=10,loc='left',y=1.04)
   except AttributeError:
     # Python 2 failing with attribute "loc".
     ax5.set_title("Antenna patterns: in use (lines) and AIS-measured (circles)", fontsize=10,horizontalalignment='left',y=1.04)

   # ...Plot the AIS patterns.

   # ......Interpolate to finer spacing if requested (in the case of gappier
   # AIS patterns, interpolation can be distracting, as large gaps get filled).
   doInterp = cfgParser.getboolean('plots','doInterp')
   if doInterp:
       interpIncrement_degrees = cfgParser.getfloat('plots','interpIncrement_degrees')
       interpTo_degrees = np.arange(0,360,interpIncrement_degrees)
       interpTo_radians = np.array([radians(x) for x in interpTo_degrees])

       isNotNan = ~np.isnan(smoothedAmp1)
       realY = np.interp(interpTo_degrees, np.real(smoothedBearing1[isNotNan]), smoothedAmp1[isNotNan])
       imagY = np.interp(interpTo_degrees, np.imag(smoothedBearing1[isNotNan]), smoothedAmp1[isNotNan])

       interpAmp1 = []
       for x in map(complex,realY,imagY):
           interpAmp1.append(x)
       interpAmp1 = np.array(interpAmp1)

       isNotNan = ~np.isnan(smoothedAmp2)
       realY = np.interp(interpTo_degrees, np.real(smoothedBearing2[isNotNan]), smoothedAmp2[isNotNan])
       imagY = np.interp(interpTo_degrees, np.imag(smoothedBearing2[isNotNan]), smoothedAmp2[isNotNan])

       interpAmp2 = []
       for x in map(complex,realY,imagY):
           interpAmp2.append(x)
       interpAmp2 = np.array(interpAmp2)

       aisRadians = interpTo_radians
       aisAmp1 = interpAmp1
       aisAmp2 = interpAmp2

       #ax5.scatter(x=interpTo_radians,y=np.abs(interpAmp1),marker='o',s=20,edgecolors=LOOP1_COLOUR, facecolors='none')
       #ax5.scatter(x=interpTo_radians,y=np.abs(interpAmp2),marker='o',s=20,edgecolors=LOOP2_COLOUR, facecolors='none')
   else:
       smoothedBearing_radians1 = np.array([radians(x) for x in smoothedBearing1])
       smoothedBearing_radians2 = np.array([radians(x) for x in smoothedBearing2])
       #ax5.scatter(x=smoothedBearing_radians,y=np.abs(smoothedAmp1),marker='o',s=20,edgecolors=LOOP1_COLOUR, facecolors='none')
       #ax5.scatter(x=smoothedBearing_radians,y=np.abs(smoothedAmp1),marker='o',s=20,edgecolors=LOOP2_COLOUR, facecolors='none')
       aisRadians1 = smoothedBearing_radians1
       aisRadians2 = smoothedBearing_radians2
       aisAmp1 = smoothedAmp1
       aisAmp2 = smoothedAmp2

   ax5.scatter(x=aisRadians1,y=np.abs(aisAmp1),marker='o',s=20,edgecolors=LOOP1_COLOUR, facecolors='none')
   ax5.scatter(x=aisRadians2,y=np.abs(aisAmp2),marker='o',s=20,edgecolors=LOOP2_COLOUR, facecolors='none')

   # ...Plot the pattern currently in use.

   # ......Older versions of numpy don't have nanmedian, nanmean, etc. For backwards
   # compatibility, call a de-nanning function before finding median or mean.

   # ......Scale the pattern in use to be similar amplitude to AIS pattern. Use
   # a single scale factor for both loops.
   #loop1Pattern_median = np.nanmedian(np.absolute(measAmp1))
   #loop2Pattern_median = np.nanmedian(np.absolute(measAmp2))

   denannedX = denan(np.absolute(measAmp1))
   if len(denannedX) > 0:
     loop1Pattern_median = np.median(denannedX)
   else:
     # Scaling factor of 1.
     loop1Pattern_median = 1

   denannedX = denan(np.absolute(measAmp2))
   if len(denannedX) > 0:
     loop2Pattern_median = np.median(np.absolute(denannedX))
   else:
     # Scaling factor of 1.
     loop2Pattern_median = 1

   if np.all(np.isnan(np.absolute(smoothedAmp1))):
       # Can't find median if all NaN.
       A13_median = loop1Pattern_median
   else:
       #A13_median = np.nanmedian(np.absolute(smoothedAmp1))
       denannedX = denan(np.absolute(smoothedAmp1))
       if len(denannedX) > 0:
         A13_median = np.median(np.absolute(denannedX))
       else:
         A13_median = 1

   if np.all(np.isnan(np.absolute(smoothedAmp2))):
       # Can't find median if all NaN.
       A23_median = loop2Pattern_median
   else:
       #A23_median = np.nanmedian(np.absolute(smoothedAmp2))
       denannedX = denan(np.absolute(smoothedAmp2))
       if len(denannedX) > 0:
         A23_median = np.median(np.absolute(denannedX))
       else:
         A23_median = 1

   #scaleFac = np.nanmean(np.array([A13_median/loop1Pattern_median, A23_median/loop2Pattern_median]))
   scaleFac = np.mean(np.array([A13_median/loop1Pattern_median, A23_median/loop2Pattern_median]))
   scaledLoop1Pattern = np.absolute(measAmp1) * scaleFac
   scaledLoop2Pattern = np.absolute(measAmp2) * scaleFac

   measBearing_radians1 = np.array([radians(x) for x in measBearing1])
   measBearing_radians2 = np.array([radians(x) for x in measBearing2])

   # Use plot, rather than scatter (assumes that the pattern in use has been
   # nicely smoothed and interpolated to fine increments, so drawing lines
   # rather than independent points won't result in a confusing tangle.
   #ax5.scatter(x=measBearing_radians, y=scaledLoop1Pattern, marker='.',s=2,color=LOOP1_COLOUR)
   #ax5.scatter(x=measBearing_radians, y=scaledLoop2Pattern, marker='.',s=2,color=LOOP2_COLOUR)
   ax5.plot(measBearing_radians1, scaledLoop1Pattern, linestyle='-', marker=None, color=LOOP1_COLOUR)
   ax5.plot(measBearing_radians2, scaledLoop2Pattern, linestyle='-', marker=None, color=LOOP2_COLOUR)

   if np.isnan(numDaysAveraged):
       # Called with NaN for this value if a single-day plot.
       titleStr = "Single-day " + stationName.upper() + " AIS pattern, " + str(date)
   else:
       titleStr = stationName.upper() + " AIS pattern, " + str(date) + " back-averaged " + str(numDaysAveraged) + " days (" + cfgParser.get('plots','timeAveragePlotLength_days') + " requested)"
   #ax5.set_title(titleStr, fontsize=12,y=1.07)
   fig.suptitle(titleStr, fontsize=14)

   gs.update(wspace=0.3, hspace=0.2)
   return fig

################################################################
def generate_plots(cfgParser,loopObj,thisDay,loop1PatternHistory,loop2PatternHistory):
     # Calls the plotting routines.
     figureOutputDir = cfgParser.get('files','figureOutputDir')
     if not os.path.exists(figureOutputDir):
       try:
          os.makedirs(figureOutputDir)
       except:
          raise Exception("Failed to create figure output directory.")

     patternFile = cfgParser.get('files','pattFile')
     loop1MeasPattern, loop2MeasPattern, measMetaData = read_pattern_file(patternFile)
     print("Pattern file", patternFile, "read in.")

     figHndls = []
     outFileNames = []

     # Plot the raw and filtered data to demonstrate the effects of the
     # filters on this day's data.
     fig1 = plot_filter_results(loopObj,thisDay,cfgParser)
     figName = "plot_filter_results_" + stationName.upper() + "_" + str(thisDay) + ".png"
     figName = os.path.join(figureOutputDir,figName)
     figHndls.append(fig1)
     outFileNames.append(figName)

     # Plot the day's AIS pattern along with the boat pattern.
     # ...Create smoothed versions of the patterns.
     smoothedBinSize_degrees = cfgParser.getfloat('processing','smoothedBinSize_degrees')
     smoothedA13pattern = Pattern(loopObj.get_filtered('targetBearing'),loopObj.get_filtered('A13'),smoothedBinSize_degrees)
     smoothedA23pattern = Pattern(loopObj.get_filtered('targetBearing'),loopObj.get_filtered('A23'),smoothedBinSize_degrees)

     # ...a NaN value here denotes a single day of data (no averaging). The
     # effect is the same as using a value of 1, but the axis title is
     # different.
     numDaysAveraged = np.nan
     #fig2 = plot_AIS_pattern(A13pattern,A23pattern,smoothedA13pattern,smoothedA23pattern, loop1MeasPattern,loop2MeasPattern,thisDay,numDaysAveraged,cfgParser)
     fig2 = plot_AIS_pattern(loop1PatternHistory.patterns[thisDay],loop2PatternHistory.patterns[thisDay],smoothedA13pattern,smoothedA23pattern, loop1MeasPattern,loop2MeasPattern,thisDay,numDaysAveraged,cfgParser)
     figName = "plot_AIS_pattern_" + stationName.upper() + "_" + str(thisDay) + ".png"
     figName = os.path.join(figureOutputDir,figName)
     figHndls.append(fig2)
     outFileNames.append(figName)

     # Plot a time-averaged pattern.
     timeAveragePlotLength_days = cfgParser.getint('plots','timeAveragePlotLength_days')
     loop1PatternHistorySubset = loop1PatternHistory.subset(timeAveragePlotLength_days,lastDate)
     loop2PatternHistorySubset = loop2PatternHistory.subset(timeAveragePlotLength_days,lastDate)
     fig3 = plot_AIS_history(loop1PatternHistorySubset,loop2PatternHistorySubset,loop1MeasPattern,loop2MeasPattern,lastDate,timeAveragePlotLength_days,cfgParser)
     figName = "plot_AIS_pattern" + "_" + stationName.upper() + "_history_" + str(timeAveragePlotLength_days) + "_" + str(thisDay) + ".png"
     figName = os.path.join(figureOutputDir,figName)
     figHndls.append(fig3)
     outFileNames.append(figName)

     # Print or display the figures generated from this loop file.

     # ...Must be no "show" or "draw" command if run as headless script.
     headless = cfgParser.getboolean('plots','headless')

     if headless:
         print(scriptName, "--Running in headless mode.")
         for figHndl,outFileName in zip(figHndls,outFileNames):
             print("Saving figure to",outFileName)
             figHndl.savefig(outFileName,dpi=200)
             plt.close()
     else:
         print(scriptName, "--Running in NON-headless mode.")
         plt.show()

################################################################
def get_config(cfgFile):
  cfgParser = configparser.ConfigParser()
  cfgParser.optionxform=str
  try:
    found = cfgParser.read(os.path.expanduser(cfgFile))
    if len(found) < 1:
      raise DailyAisError(".ini file not found")
    return cfgParser
  except Exception as e:
    print(e)
    raise DailyAisError("Failed to read .ini file", cfgFile)

################################################################
def get_filters(cfgParser):
    filters = {}
    filters['minSNR'] = cfgParser.getfloat('filters','minSNR')
    filters['maxSNR'] = cfgParser.getfloat('filters','maxSNR')
    filters['minRangeWidth'] = cfgParser.getfloat('filters','minRangeWidth')
    filters['maxRangeWidth'] = cfgParser.getfloat('filters','maxRangeWidth')
    filters['minDoppWidth'] = cfgParser.getfloat('filters','minDoppWidth')
    filters['maxDoppWidth'] = cfgParser.getfloat('filters','maxDoppWidth')
    filters['minRange'] = cfgParser.getfloat('filters','minRange')
    filters['maxRange'] = cfgParser.getfloat('filters','maxRange')
    filters['minVesselSpeed'] = cfgParser.getfloat('filters','minVesselSpeed')
    filters['maxVesselSpeed'] = cfgParser.getfloat('filters','maxVesselSpeed')
    filters['nearBragg'] = cfgParser.getfloat('filters','nearBragg')

    filters['bearingStart'] = []
    for x in cfgParser.get('filters','bearingStart').split():
        filters['bearingStart'].append(float(x))

    filters['bearingEnd'] = []
    for x in cfgParser.get('filters','bearingEnd').split():
        filters['bearingEnd'].append(float(x))
    return filters

################################################################
def parse_loop_name(filename):
    filePath = os.path.dirname(filename)
    baseName = os.path.basename(filename)
    baseBaseName, fileExtension = os.path.splitext(baseName)
    strs = baseBaseName.split('_')
    fileStationName = strs[1]
    dateStr = "-".join([strs[2],strs[3],strs[4],strs[5]])
    fileDatetime = dt.datetime.strptime(dateStr,'%Y-%m-%d-%M%S')
    return filePath, fileStationName, fileDatetime

################################################################
def denan(x):
  # For backwards compatibility with older versions of numpy,
  # provide this function to remove NaNs from arrays before
  # finding median, mean, etc. (older numpy has no nanmedian,
  # nanmean, etc.). This function is much less flexible and
  # probably less robust than the official numpy version,
  # but it should do for this application.
  denannedX = x[~np.isnan(x)]
  return denannedX

################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store_true', help="Force re-read of all loop files.")
    args, whatsLeft = parser.parse_known_args()

    if len(whatsLeft) > 0:
       cfgFile = whatsLeft[0]
    else:
      raise DailyAisError("Must specify a config file")

    if args.f:
      doForce = True
      print("doForce set to True. Will read in all available loop files.")
    else:
      print("doForce set to False. Will read in only those loop files that are not already represented in the pattern history file.")
      doForce = False

    cfgParser = get_config(cfgFile);
    DOMAP = cfgParser.getboolean('map','DOMAP')

    if DOMAP:
        from mpl_toolkits.basemap import Basemap

    today = dt.datetime.combine(dt.date.today(), dt.datetime.min.time())

    loopDir = cfgParser.get('files','loopDir')
    dataOutputDir = cfgParser.get('files','dataOutputDir')

    if not os.path.exists(dataOutputDir):
        try:
          os.makedirs(dataOutputDir)
        except:
          raise Exception("Failed to create data output directory.")

    stationName = cfgParser.get('station','name')
    binSize_degrees = cfgParser.getfloat('processing','binSize_degrees')

    # Load the bearing-averaged pattern history, if it exists. Do not 
    # test for existence of file directly--shelve sometimes adds ".db"
    # or ".dat" to the end of filenames.
    patternHistoryFile = cfgParser.get('files','patternHistoryFile')
    try:
      shelf = shelve.open(patternHistoryFile, flag='r')
      print("Successfully opened patternHistoryFile.")
    except ShelfOpenError:
      print("Failed to open patternHistoryFile. Create an empty one...")

      # History does not exist. Create an empty one.
      # 2018-10-10, kpb--Mysterious overwriting of existing shelf file occurred. Don't really understand
      # why, but will try to prevent it happening again by testing for existence of any file with a
      # matching name.
      if len(glob.glob(patternHistoryFile + '*')) > 0:
        raise DailyAisError("patternHistoryFile already exists. Cannot create a new one until it is deleted.")

      try:
        shelf = shelve.open(patternHistoryFile, flag='c')
        loop1PatternHistory = PatternHistory()
        loop2PatternHistory = PatternHistory()
        shelf['loop1PatternHistory'] = loop1PatternHistory
        shelf['loop2PatternHistory'] = loop2PatternHistory
      except:
        raise DailyAisError("Failed to create patternHistoryFile.")
      
    # 2020-10-09, kpb. Workaround for apparent python bug that results in KeyError, even 
    # when "'loop2PatternHistory' in shelf" returns True. Listing the keys is a fast way
    # to cause an exception, detecting a corrupted shelve object and bail out gracefully.
    try:
       tmp = list(shelf.keys())
    except:
       shelf.close()
       raise DailyAisError("Failed to read keys from shelved patternHistoryFile. File may be corrupt. Delete it and try again. ")

    try:
      loop1PatternHistory = shelf['loop1PatternHistory']
      loop2PatternHistory = shelf['loop2PatternHistory']
      shelf.close()
    except:
      raise DailyAisError("Failed to load shelved patternHistoryFile")

    # Make list of loop files that are not already in the history.
    loopFilesToProcess = []
    searchStr = os.path.join(loopDir,'*.loop')
    allLoopFiles = glob.glob(os.path.join(loopDir,'*.loop'))

    for loopFile in allLoopFiles:
        filePath, fileStationName, fileDatetime = parse_loop_name(loopFile)
        if fileDatetime < today:
           # The nominal date of this loop file is from before today, so it can safely
           # be processed (processing a loop file from TODAY could conceivably result
           # in a partial day's worth of loop data being read in; the day would be 
           # flagged as having already been processed and the remainder of the day's
           # loop data would never get read).
           if doForce or fileDatetime.date() not in loop1PatternHistory.patterns.keys():
               # Either this loop file has not already been processed OR the operator has instructed
               # that all loop files are to be processed regardless. Either way, add this file to the
               # list of loop files to process.
               loopFilesToProcess.append(loopFile)
    loopFilesToProcess.sort()
    numLoops = len(loopFilesToProcess)
    print(len(allLoopFiles),"loop files available. Will process",numLoops, "of them.")

    # Process the daily loop files.
    loopCount = 1
    historyModified = False

    for loopFile in loopFilesToProcess:
        print("Processing loop file",loopCount,"of",numLoops,":",loopFile)
        loopCount = loopCount + 1

        # Read in this day's loop file.
        try:
            loopObj = Codar_loop(loopFile)
            thisDay = loopObj.calc_time()[0].date()
        except:
            print("Failed to read in loop file",loopFile,". Skipping...")
            continue

        loopTime = loopObj.time[0]
        lastDate = loopTime[0].date()

        # Filter this day's loop data according to parameters in Chad Whelan's
        # Matlab routine aispatt2.m.
        filters = get_filters(cfgParser)
        loopObj.set_filters(filters)
        loopObj.apply_filters()

        # Create patterns from this day's filtered loop data.
        A13pattern = Pattern(loopObj.get_filtered('targetBearing'),loopObj.get_filtered('A13'),binSize_degrees)
        A23pattern = Pattern(loopObj.get_filtered('targetBearing'),loopObj.get_filtered('A23'),binSize_degrees)

        # Add this day's bearing averaged patterns (loop 1 and loop 2) to the history.
        loop1PatternHistory.add(A13pattern,thisDay)
        loop2PatternHistory.add(A23pattern,thisDay)
        print("thisDay",thisDay,"pattern added to pattern history")
        historyModified = True

        # Generate the plots.
        generate_plots(cfgParser,loopObj,thisDay,loop1PatternHistory,loop2PatternHistory)

        # Check to see if it is time to store the bearing-averaged pattern
        # history in a yearly file.
        years = np.unique([x.year for x in loop1PatternHistory.patterns.keys()])
        for year in years:
            if year < np.max(years):
                # Can't simply check for existence of shelf file because shelve sometimes adds ".db" or ".dat" to the end
                # of the user-specified filename. Do a glob instead. n.b., shelve may create more than one shelf file with the same basename.
                #thisYearHistoryFile = os.path.join(dataOutputDir,"patternHistory_" + stationName.upper() + "_" + str(year) + ".shv")
                thisYearHistoryFiles = glob.glob(os.path.join(dataOutputDir,"patternHistory_" + stationName.upper() + "_" + str(year) + "*"))
                if len(thisYearHistoryFiles) < 1:
                    # No shelf file matching the year exists. Create one.
                    thisYearHistoryFile = os.path.join(dataOutputDir,"patternHistory_" + stationName.upper() + "_" + str(year) + ".shv")
                    loop1YearSubset = loop1PatternHistory.extract_year(year)
                    loop2YearSubset = loop2PatternHistory.extract_year(year)
                    print("Saving the pattern history for year",year,"...")
                    try:
                        shelf = shelve.open(thisYearHistoryFile, flag='c')
                        shelf['loop1PatternHistory'] = loop1YearSubset
                        shelf['loop2PatternHistory'] = loop2YearSubset
                        shelf.close()
                    except:
                        raise DailyAisError("Failed to save history file for year",year)

                    print("Saving yearly pattern history file complete.")

        # Retain only a set number of days worth of pattern history in the current
        # history file. This number of days should be greater than a year, so that
        # the data will not have been deleted when the yearly history file is
        # made.
        DAYS_TO_RETAIN = 370
        #print("Pre trim",loop1PatternHistory.patterns.keys())
        loop1PatternHistory.trim(DAYS_TO_RETAIN,lastDate)
        loop2PatternHistory.trim(DAYS_TO_RETAIN,lastDate)
        #print("Post trim",loop1PatternHistory.patterns.keys())

    # Save the newly-modified pattern history if any new loop files were read in.
    if historyModified:
        print("Saving the bearing-averaged pattern history...")
        try:
            parentDir = os.path.dirname(patternHistoryFile)
            if not os.path.exists(parentDir):
                os.makedirs(parentDir)

            shelf = shelve.open(patternHistoryFile, flag='w')
            print("open successful")
            shelf['loop1PatternHistory'] = loop1PatternHistory
            shelf['loop2PatternHistory'] = loop2PatternHistory
            shelf.close()
        except:
            raise  DailyAisError("Failed to save modified patternHistoryFile.")

        print("Saving of bearing-averaged history complete. File",patternHistoryFile)

