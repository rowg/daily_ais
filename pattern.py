''' 
REVISION HISTORY

2018-09-12, kpb--Backwards compatibility for older numpy.

2018-10-05, kpb--Test to prevent finding average of non-finite values.

'''

import numpy as np
import matplotlib.pyplot as plt
from math import radians
import datetime as dt
import string
import copy
import shelve

class PatternError(Exception):
   """Basic exception for Pattern class"""

################################################################
class Pattern(object):
    # A Pattern object represents an antenna pattern containing complex amplitudes
    # binned by bearing. Each binned complex amplitude is associated with a weight, which
    # indicates the number of data points that contributed to the binned amplitude.
    def __init__(self,bearings,amplitudes,binSizeDeg,weights=None):
        self.binSizeDeg = binSizeDeg
        #bearingBins = np.arange(binSizeDeg,360+binSizeDeg,binSizeDeg)
        bearingBins = np.arange(0,360,binSizeDeg)
        if weights is None:
            weights = np.array([1]*len(bearings))
        averages,numPts = self.average_by_bearing_bins(amplitudes,bearings,bearingBins,weights)
        self.amplitudes = averages
        self.weights = numPts

    def get_arrays(self):
        # Returns bearings, amplitudes and weights of the pattern, excluding
        # those with weights of zero (i.e., for bearings with no data points).
        bearings = []
        amplitudes = []
        weights = []
        for thisBearing in self.amplitudes.keys():
            if self.weights[thisBearing] != 0:
                bearings.append(thisBearing)
                amplitudes.append(self.amplitudes[thisBearing])
                weights.append(self.weights[thisBearing])
        bearings = np.array(bearings)
        amplitudes = np.array(amplitudes)
        weights = np.array(weights)
        return bearings, amplitudes, weights
        
    def average_by_bearing_bins(self,x,bearings,bearingBins,weights=None):
        # Bins array "x" by array "bearings" and finds average for each bin. Bin values given
        # by array "bearings" represent the CENTRE of each bin.
        #
        # If optional input argument weights (numpy array of integers, the same
        # length as x and bearings) is specified, then the average is weighted by the
        # specified value; this prevents an average of averages being unduly
        # affected by an outlier caused by a single data point, for example.
        if weights is None:
            print("WEIGHTS IS NONE")
            weights = np.array([1]*len(x))
        if not (np.size(x) == np.size(bearings)) and (np.size(x) == np.size(weights)):
            raise PatternError("Arrays x, bearings and weights must be the same size.")

        # Numpy's digitize() bin values represent the UPPER limit of each bin. Shift
        # the bins by half a bin width so that index returned is to the CENTRE of
        # each bin instead.
        if len(np.unique(np.diff(bearingBins))) != 1:
            raise PatternError("Bearing bins should be of uniform width.")
        binWidth = np.unique(np.diff(bearingBins))[0]
        shiftedBins = bearingBins + 0.5*binWidth

        # ...Binning bearings near the 359 degree to 0 degree wrap-around is a
        # bit tricky. Example: assume a bin width of 1 degree. We would want
        # a bearing value of 359.8 to be binned in the bin centred at 0 degrees, 
        # because 359.8 is closer to 0 degrees than it is to the last bin,
        # which is centred at 359 degrees. Achieve this by wrapping bearings
        # of greater than 359.5 degrees to negative values.
        shiftedBearings = np.copy(bearings) # Have to do this or python will change the original when bearings is wrapped.
        shiftedBearings[shiftedBearings>shiftedBins[-1]] = shiftedBearings[shiftedBearings>shiftedBins[-1]] - 360
        binInd = np.digitize(shiftedBearings, shiftedBins)
        averages = dict(zip(bearingBins,[complex(np.nan,np.nan)]*len(bearingBins)))
        numPts = dict(zip(bearingBins,[0]*len(bearingBins)))
        for i in range(len(bearingBins)):
          thisBinVal = bearingBins[i]
          thisBinInds = binInd==i
          numThisBin = np.size(np.nonzero(thisBinInds))
          if numThisBin>0:
             thisBinX = x[thisBinInds]
             thisBinWeights = weights[thisBinInds]
             #print("thisBinVal is",thisBinVal)
             #print("thisBinX is",thisBinX,"magnitude",abs(thisBinX))
             #print("thisBinWeights is",thisBinWeights)

             # 2018-10-05, kpb--LOOP_VCOL_2018_05_11_0000.loop has "inf" and "nan" in some data columns,
             # causing an error. Only try to find the average of finite values.
             isFinite = np.isfinite(abs(thisBinX))
             thisBinX = thisBinX[isFinite]
             thisBinWeights = thisBinWeights[isFinite]
             if len(thisBinX) > 0:
               thisBinWeightedVals = thisBinWeights * thisBinX
               thisBinWeightedSum = np.sum(thisBinWeightedVals)
               thisBinWeight = np.sum(thisBinWeights)
               averages[thisBinVal] = thisBinWeightedSum / thisBinWeight
               numPts[thisBinVal] = thisBinWeight
        return averages,numPts

################################################################
class PatternHistory(object):
    # A PatternHistory object contains multiple Pattern objects as
    # a dictionary keyed by dates.
    def __init__(self):
        self.patterns = {}

    def add(self,newPattern,date):
        #print("adding")
        self.patterns[date] = newPattern
        
    def get_arrays(self,dates=None):
        #print("getting")
        bearings = []
        amplitudes = []
        weights = []
        if dates is None:
            dates = self.patterns.keys()
        #print("dates is",dates)
        for date in dates:
            #print("date is",date)
            if date in self.patterns.keys():
              thisPatternBearings, thisPatternAmplitudes, thisPatternWeights = self.patterns[date].get_arrays()
              bearings.extend(thisPatternBearings[:])
              amplitudes.extend(thisPatternAmplitudes[:])
              weights.extend(thisPatternWeights[:])
        bearings = np.array(bearings)
        amplitudes = np.array(amplitudes)
        weights = np.array(weights)
        return bearings, amplitudes, weights

    def get_dates(self,numDates,lastDate):
        # Returns up to a maximum of numDates dates going backwards starting 
        # and including) lastDate. Only dates that contain patterns are returned.
        datesRequested = [lastDate - dt.timedelta(days=1) * i for i in range(numDates)]
        availableDates = np.intersect1d(datesRequested,list(self.patterns.keys()))
        #print("list(self.patterns.keys())",list(self.patterns.keys()))
        #print("availableDates",availableDates)
        #print("numDates",numDates,"lastDate",lastDate,"datesRequested",datesRequested)
        return availableDates

    def subset(self,numDates,lastDate):
        # Returns subset of the PatternHistory object without changing the
        # original object. The subset is a maximum of numDates days long; 
        historySubset = PatternHistory()
        availableDates = self.get_dates(numDates,lastDate)
        #print("subset numDates",numDates,"lastDate",lastDate,"availableDates",availableDates)
        for thisDate in availableDates:
            historySubset.add(self.patterns[thisDate],thisDate)
        return historySubset
    
    def extract_year(self,year):
        # Returns a subset of the PatternHistory containing all the patterns
        # from the specified year.
        yearSubset = PatternHistory()
        for thisDate in self.patterns.keys():
            if thisDate.year == year:
                yearSubset.add(self.patterns[thisDate],thisDate)
        return yearSubset
    
    def trim(self,numDates,lastDate):
        # Modifies the original PatternHistory object by removing patterns
        # that are not in the specified date range.
        datesToKeep = self.get_dates(numDates,lastDate)
        #print("datesToKeep",datesToKeep)
        #print("list of keys",list(self.patterns.keys()))
        # Older versions of numpy throw an exception when setdiff1d is applied
        # to datetimes. 
        try:
          datesToDelete = np.setdiff1d(np.array(list(self.patterns.keys())),datesToKeep)
        except TypeError:
          # This approach is probably slower, but should work in python 2:
          datesToDelete = []
          for key in self.patterns.keys():
            if key in datesToKeep:
               datesToDelete.append(key)
          datesToDelete = np.array(datesToDelete)

        for thisDate in datesToDelete:
           #print("Delete thisDate",thisDate)
           del self.patterns[thisDate]

################################################################
def read_pattern_file(inFile):
  # Reads in a MeasPattern.txt or IdealPattern.txt codar file.
  inFid = open(inFile,'r')
  line = 'dummy'

  # Read in first line. This gives number of bearing values.
  line = inFid.readline().strip()
  if len(line)<=0:
     print("Pattern file is empty. Aborting.")
     raise PatternError("Empty pattern file")

  numBearings = int(line)

  # Don't know why, but there are 9 data tables in MeasPattern.txt, each
  # with numBearings elements. Several of these are all-zero. Same applies
  # to IdealPattern.txt. Will try reading in and plotting to see if I can
  # figure them out.
  numTables = 9
  numCells = numTables*numBearings
  tablesCompleted = False
  data = []

  # IdealPattern.txt files don't contain all the metadata that MeasPattern.txt
  # files do. Initialise all the metadata variables to None.
  metaData = {}
  metaData['AntennaBearing'] = None
  metaData['DateTime'] = None
  metaData['AmplitudeFactor1'] = None
  metaData['AmplitudeFactor2'] = None
  metaData['SiteCode'] = None
  metaData['SiteLat'] = None
  metaData['SiteLon'] = None
  metaData['DegreeResolution'] = None
  metaData['DegreeSmoothing'] = None
  metaData['UUID'] = None
  metaData['PhaseCorrection1'] = None
  metaData['PhaseCorrection2'] = None
  
  cellCount = 0
  while len(line)>0:
     line = inFid.readline().strip() 
     if len(line)>0:
         if not tablesCompleted:
             cols = line.split()      
             for col in cols:
                 data.append(float(col))
                 cellCount = cellCount + 1
                 if cellCount >= numCells:
                     tablesCompleted = True

         else:
             cols = line.split('!')
             colName = cols[-1].strip(string.punctuation + string.whitespace)
             if colName == 'Antenna Bearing':
                 AntennaBearing = float(cols[0])
                 metaData['AntennaBearing'] = AntennaBearing
             elif colName == 'Amplitude Factors':
                 strs = cols[0].split()
                 metaData['AmplitudeFactor1'] = float(strs[0])
                 metaData['AmplitudeFactor2'] = float(strs[1])
             elif colName == 'Site Code':
                 metaData['SiteCode'] = cols[0].strip(string.punctuation + string.whitespace)
             elif colName == 'Site Lat Lon':
                 strs = cols[0].split()
                 metaData['SiteLat'] = float(strs[0])
                 metaData['SiteLon'] = float(strs[1])
             elif colName == 'Degree Resolution':
                 metaData['DegreeResolution'] = float(cols[0])
             elif colName == 'Degree Smoothing':
                 metaData['DegreeSmoothing'] = float(cols[0])
             elif colName == 'Date Year Mo Day Hr Mn Sec':
                 # Not sure how consistent date format is (extra spaces?),
                 # so convert to individual integers and THEN to datetime.
                 strs = cols[0].split()
                 metaData['DateTime'] = dt.datetime(int(strs[0]),int(strs[1]), int(strs[2]), int(strs[3]), int(strs[4]), int(strs[5]), 0)
             elif colName == 'UUID':
                 metaData['UUID'] = cols[0].strip(string.punctuation + string.whitespace)
             elif colName == 'Phase Corrections':
                 strs = cols[0].split()
                 metaData['PhaseCorrection1'] = float(strs[0])
                 metaData['PhaseCorrection2'] = float(strs[1])

  inFid.close()

  data = np.array(data)
  data = np.reshape(data,(numTables,-1))
  ccwBearings = data[0,:]

  # Trial and error show that data columns are 1 and 3 for loop 1, and 5 and 7 for loop 2.
  loop1 = []
  for x in map(complex,data[1,:],data[3,:]):
      loop1.append(x)
  loop1 = np.array(loop1)

  loop2 = []
  for x in map(complex,data[5,:],data[7,:]):
      loop2.append(x)
  loop2 = np.array(loop2)

  # Bearings in Seasonde .patt files are given as degrees clockwise from
  # the True antenna bearing. In contrast, bearings in the 
  # MeasPattern.txt and IdealPattern.txt are apparently given as degrees
  # COUNTER-clockwise from the antenna bearing (this is the only way to 
  # make my plots agree with those produced by Seasonde and to stop the
  # boat pattern appearing landward of the antenna).
  bearingsTrue = -ccwBearings + AntennaBearing
  # ...Convert negative compass directions to range 0 - 360:
  bearingsTrue[bearingsTrue<0] = bearingsTrue[bearingsTrue<0] + 360
  loop1Pattern = Pattern(bearingsTrue,loop1,1)
  loop2Pattern = Pattern(bearingsTrue,loop2,1)
  #print(loop1Pattern.amplitudes)
  return loop1Pattern, loop2Pattern, metaData
    
def plot_both_patterns(pattern1,pattern2,metaData):
    # Mostly for development/diagnostics.
    LOOP1_COLOUR = 'red'
    LOOP2_COLOUR = 'blue'

    bearings, amplitudes, weights = pattern1.get_arrays()
    bearingsRadians1 = [radians(x) for x in bearings]
    mag1 = np.absolute(amplitudes)

    bearings, amplitudes, weights = pattern2.get_arrays()
    bearingsRadians2 = [radians(x) for x in bearings]
    mag2 = np.absolute(amplitudes)

    fig = plt.figure()
    ax = plt.subplot(111, polar=True)
    ax.scatter(x=bearingsRadians1, y=mag1, marker='.',s=20,color=LOOP1_COLOUR)
    ax.scatter(x=bearingsRadians2, y=mag2, marker='.',s=20,color=LOOP2_COLOUR)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rmax(np.max(np.maximum(mag1,mag2)))

    # IdealPattern.txt files don't contain all the fields that 
    # MeasPattern.txt files do.
    if metaData['SiteCode'] is not None:
        SiteCode = metaData['SiteCode']
    else:
        SiteCode = "Ideal"

    if metaData['DateTime'] is not None:
        DateStr = metaData['DateTime'].strftime('%Y-%m-%d %H:%M:%S')
    else:
        DateStr = ""

    plt.title(SiteCode + " Antenna Pattern, " + DateStr,y=1.08, va='bottom')

if __name__ == '__main__':
  print("Testdrive starting...")

  # Create a pattern object from some dummy bearing and amplitudes.

  # ...Start with poorly-sorted numpy arrays of bearings and complex amplitudes with 
  # multiple amplitude values for a single bearing like you might get from a loop file.
  myBearings = np.array([1, 0.4, 2, 3, 4, 3, 1, 3, 245, 245, 246, 247, 359.4, 359.7])
  myAmplitudes = np.array([np.complex(1,9), np.complex(2,5), np.complex(4,5), np.complex(3,5),
                           np.complex(2,3), np.complex(2,7), np.complex(4,7), np.complex(11,35), 
                           np.complex(4,1), np.complex(1,1), np.complex(1,4), np.complex(2,3),
                           np.complex(2,5), np.complex(7,3)])
  print("")
  print("Raw data has",len(myBearings),"complex amplitudes at the following bearings:")
  print(myBearings)
  
  initBinSizeDeg = 1 
  initialPattern = Pattern(myBearings,myAmplitudes,initBinSizeDeg)
  print("")
  print("Initial pattern (n.b., the bearing value is the CENTRE of the bin):")
  for thisBearing in initialPattern.amplitudes.keys():
      if not np.isnan(initialPattern.amplitudes[thisBearing]):
          #print(thisBearing,initialPattern.amplitudes[thisBearing],initialPattern.weights[thisBearing])
          binLowerBound = thisBearing - initBinSizeDeg/2.
          binUpperBound = thisBearing + initBinSizeDeg/2.
          binStr = "Bin %005.1f (%+006.1f - %+006.1f degrees)" % (thisBearing,binLowerBound,binUpperBound)
          print(binStr,initialPattern.amplitudes[thisBearing],initialPattern.weights[thisBearing])
  print("")

  # Take the initial pattern and use it as the inputs to create a smoothed pattern.
  initialPatternBearings, initialPatternAmplitudes, initialPatternWeights = initialPattern.get_arrays()
  smoothedBinSizeDeg = 90 
  smoothedPattern = Pattern(initialPatternBearings,initialPatternAmplitudes,smoothedBinSizeDeg,initialPatternWeights)
  print("Smoothed pattern:")
  sortedBearings = list(smoothedPattern.amplitudes.keys())
  sortedBearings.sort()
    
  for thisBearing in sortedBearings:
      if not np.isnan(smoothedPattern.amplitudes[thisBearing]):
          #print(thisBearing,smoothedPattern.amplitudes[thisBearing],smoothedPattern.weights[thisBearing])
          binLowerBound = thisBearing - smoothedBinSizeDeg/2.
          binUpperBound = thisBearing + smoothedBinSizeDeg/2.
          binStr = "Bin %005.1f (%+006.1f - %+006.1f degrees)" % (thisBearing,binLowerBound,binUpperBound)
          print(binStr,smoothedPattern.amplitudes[thisBearing],smoothedPattern.weights[thisBearing])
  print("")

  # Create a pattern history.
  now = dt.datetime.now()
  loop1PatternHistory = PatternHistory()
  #loop1PatternHistory.add(loop1MeasPattern,now.date())

  # Add multiple dummy patterns, each different, to the history.
  dayCount = 0
  historyDays = []
  for i in [1,10,100,1000]:
    thisDateTime = now + dt.timedelta(days=dayCount)
    historyDays.append(thisDateTime.date())
    dayCount = dayCount+1
    newPattern = copy.deepcopy(initialPattern) 
    for thisKey in newPattern.amplitudes.keys():
      newPattern.amplitudes[thisKey] =  i*newPattern.amplitudes[thisKey]
    loop1PatternHistory.add(newPattern,thisDateTime.date())
  historyDays = np.array(historyDays)
        
  # Extract last three of the four days of data from the pattern history and 
  # use the resulting numpy arrays as inputs to create a new, time-averaged,
  # bearing-binned pattern.
  dates = historyDays[1:]
  #print("dates is",dates)
  #print("historyDays is",historyDays)
  historyPatternBearings, historyPatternAmplitudes, historyPatternWeights = loop1PatternHistory.get_arrays(dates)
  #print(historyPatternBearings)
  #print("Length of history variables:",np.size(historyPatternBearings),np.size(historyPatternAmplitudes),np.size(historyPatternWeights))

  historyBinSizeDeg = 90 
  historyPattern = Pattern(historyPatternBearings,historyPatternAmplitudes,historyBinSizeDeg,historyPatternWeights)

  print("Cumulative pattern from",len(dates),"days of history:")
  sortedBearings = list(historyPattern.amplitudes.keys())
  sortedBearings.sort()
    
  for thisBearing in sortedBearings:
      if not np.isnan(historyPattern.amplitudes[thisBearing]):
          binLowerBound = thisBearing - historyBinSizeDeg/2.
          binUpperBound = thisBearing + historyBinSizeDeg/2.
          binStr = "Bin %005.1f (%+006.1f - %+006.1f degrees)" % (thisBearing,binLowerBound,binUpperBound)
          print(binStr,historyPattern.amplitudes[thisBearing],historyPattern.weights[thisBearing])
  print("")

  print("Read in a IdealPattern.txt file and plot the results:")
  inFile = "/Users/kpb/home/instruments/codarNotDropbox/AIS/automatedAIS/sampleData/vrid/IdealPattern_vrid.txt"
  loop1IdealPattern, loop2IdealPattern, idealMetaData = read_pattern_file(inFile)
  plot_both_patterns(loop1IdealPattern,loop2IdealPattern,idealMetaData)
  
  print("Read in a MeasPattern.txt file and plot the results:")
  inFile = "/Users/kpb/home/instruments/codarNotDropbox/AIS/automatedAIS/sampleData/vrid/MeasPattern_vrid.txt"
  loop1MeasPattern, loop2MeasPattern, measMetaData = read_pattern_file(inFile)
  plot_both_patterns(loop1MeasPattern,loop2MeasPattern,measMetaData)
  print("")

  print("Read in a pattern history 'shelf' file")
  patternHistoryFile = "/Users/kpb/home/instruments/codarNotDropbox/AIS/automatedAIS/sampleData/vrid/bearingAveragedPatternHistory_vrid.shv"
  
  try:
      shelf = shelve.open(patternHistoryFile, flag='r')
      loop1PatternHistory = shelf['loop1PatternHistory']
      loop2PatternHistory = shelf['loop2PatternHistory']
      shelf.close()
  except:
      raise PatternError("Unable to get 'shelved' data from specified patternHistoryFile ", patternHistoryFile)

  print("Pattern history originally of length",len(loop1PatternHistory.patterns.keys()))
  
  numDates = 5
  lastDate = dt.date(2017,12,9)
  
  loop1PatternHistorySubset = loop1PatternHistory.subset(numDates,lastDate)
  print("Pattern history subset of",len(loop1PatternHistorySubset.patterns.keys()),"days returned.")
  print("")

  print("Trim the pattern history to a length of only two dates")
  numDates = 2
  loop1PatternHistory.trim(numDates,lastDate)
  print("Pattern history now of length",len(loop1PatternHistory.patterns.keys()))
  #raise PatternError("Eject! Eject! Eject!")
  plt.show() 
  print("")
  print("Testdrive complete")
  
