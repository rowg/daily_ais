''' Codar_loop.py--Class definition for reading in SeaSonde "Loop" files 
containing Antenna Pattern Measurement (APM) data.

REVISION HISTORY

2018-10-02, kpb--Major bugfix. Was using numpy's complex() to generate
complex amplitudes from magnitude and phase, but complex() in fact
takes Real and Imaginary components as inputs.

2018-10-05, kpb--Removed unnecessary call to read_ctf() (already called by CodarTableFormat's __init__()).

'''

from CodarTableFormat import CodarTableFormat, CodarTableFormatError
import datetime as dt
import numpy as np
from math import radians

class CodarLoopError(CodarTableFormatError):
   """Basic exception for Codar_loop class"""

class Codar_loop(CodarTableFormat):

    ################################################################
    # Override the parent class' __init__().
    def __init__(self,inFile,userDefinedColNames=None):
        # The child's __init__() calls the parent's __init__()
        # using super():
        super(Codar_loop,self).__init__(inFile,userDefinedColNames)
        #self.read_ctf()
        self.time = {}
        self.time[0] = self.calc_time()
        
        # Loop files can contain target bearings of 360., sometimes with a little rounding error in addition. Correct this to zeroes.
        tableNum = 0
        bearing = np.array(self.get_column('TRGB',tableNum))
        bearing[bearing>=360] = 0
        bearing = list(bearing)
        self.set_column('TRGB',bearing,tableNum)
        self.TransmitCenterFreqMHz = float(self.data['fileHeader']['variables']['TransmitCenterFreqMHz'])
        self.loop2pattern_variables()
        self.filters = {}
        self.filters['minSNR'] = -np.Inf
        self.filters['maxSNR'] = np.Inf
        self.filters['minRangeWidth'] = -np.Inf
        self.filters['maxRangeWidth'] = np.Inf 
        self.filters['minDoppWidth'] = -np.Inf
        self.filters['maxDoppWidth'] = np.Inf
        self.filters['minRange'] = -np.Inf
        self.filters['maxRange'] = np.Inf
        self.filters['minVesselSpeed'] = -np.Inf
        self.filters['maxVesselSpeed'] = np.Inf
        self.filters['nearBragg'] = 0
        self.filters['bearingStart'] = [0]
        self.filters['bearingEnd'] = [360]

    ################################################################
    # Override the parent class' calc_time() to match the time format in Codar loop files.
    def calc_time(self, tableNum=None):
      # Ignore the tableNum argument; there is only one table in .loop files.
      tableNum = 0

      # Time format is space-delimited
      #   YYYYMMDD HHMMSS
      # with column names
      #   DATE TIME
      dateFlt = self.get_column('DATE',tableNum)
      timeFlt = self.get_column('TIME',tableNum)
      ctfTime = []
      for i in range(len(dateFlt)):
        dateTimeStr = "%d%06d" % (dateFlt[i],timeFlt[i])
        ctfTime.append(dt.datetime.strptime(dateTimeStr,'%Y%m%d%H%M%S'))
      return ctfTime

    ################################################################
    def get_filtered(self,fieldName):
        #return self.A13[self.doKeep['master']], self.A23[self.doKeep['master']], self.targetBearing[self.doKeep['master']]
        # getattr(self.A13,self.doKeep['master']
        return getattr(self,fieldName)[self.doKeep['master']]

    ################################################################
    def set_filters(self,filters):
        filterNames = ['minSNR', 'maxSNR','minRange','maxRange','minVesselSpeed','maxVesselSpeed','minRangeWidth','maxRangeWidth', 'minDoppWidth','maxDoppWidth','nearBragg','bearingStart', 'bearingEnd']
        for thisFilter in filters.keys():
            if thisFilter in filterNames:
                self.filters[thisFilter] = filters[thisFilter]
            else:
                raise CodarLoopError("Filter name",thisFilter,"not recognised.")

    ################################################################
    def polar2rectangular(self,mag, phase):
      return mag * np.exp(1j*phase)

    ################################################################
    def loop2pattern_variables(self):
        # Generates the complex amplitudes and bearing arrays needed to create 
        # an antenna-pattern object.
        targetBearing = np.array(self.get_column('TRGB'))
        A13magnitude = np.array(self.get_column('A13M'))
        A13phase = np.array([radians(x) for x in self.get_column('A13P')])
        A23magnitude = np.array(self.get_column('A23M'))
        A23phase = np.array([radians(x) for x in self.get_column('A23P')])

        # 2018-10-02--Major bugfix. Use new function polar2rectangular() to
        # generate complex amplitudes from magnitude and phase (numpy's
        # complex() function takes Real and Imaginary components as inputs,
        # not magnitude and phase).
        A13 = []    
        for thisComplex in map(self.polar2rectangular,A13magnitude,A13phase):
            A13.append(thisComplex)
        A13 = np.array(A13)
            
        A23 = []    
        for thisComplex in map(self.polar2rectangular,A23magnitude,A23phase):
            A23.append(thisComplex)
        A23 = np.array(A23)

        #self.data['fileHeader'] = {}
        self.A13 = A13
        self.A23 = A23
        self.targetBearing = targetBearing

    ################################################################
    def apply_filters(self):
      # If a filter value is NaN, then skip that filter.
      doKeep = {}
      
      # Chad's program looks for A3SN and either of A1SN or A2SN within limits.
      y1 = np.array(self.get_column('A1SN'))
      y2 = np.array(self.get_column('A2SN'))
      y3 = np.array(self.get_column('A3SN'))
    
      thisFiltName = 'minSNR'
      doKeep[thisFiltName] = np.array([True]*len(y1))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = np.logical_and(y3>=self.filters[thisFiltName],np.logical_or(y1>=self.filters[thisFiltName],y2>=self.filters[thisFiltName]))
          
      thisFiltName = 'maxSNR'
      doKeep[thisFiltName] = np.array([True]*len(y1))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = np.logical_and(y3<=self.filters[thisFiltName],np.logical_or(y1<=self.filters[thisFiltName],y2<=self.filters[thisFiltName]))
          
      y = np.array(self.get_column('TRGD')) # Target distance [m]
      thisFiltName = 'minRange'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y>=self.filters[thisFiltName]
    
      thisFiltName = 'maxRange'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y<=self.filters[thisFiltName]
    
      y = np.absolute(np.array(self.get_column('TRGV'))) # Target speed [m/s]
      thisFiltName = 'minVesselSpeed'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y>=self.filters[thisFiltName]
    
      thisFiltName = 'maxVesselSpeed'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y<=self.filters[thisFiltName]
    
      # Chad's Matlab: % Filter out solutions with peaks narrow in range
      #   loop(:,25) - loop(:,24) <= minRangeWidth
      # i.e., "PkWid Rng2" - "PkWid Rng1"
      # i.e.,  PKR2 - PKR1 
      y = np.array(self.get_column('PKR2')) - np.array(self.get_column('PKR1'))
    
      thisFiltName = 'minRangeWidth'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y>=self.filters[thisFiltName]
      
      thisFiltName = 'maxRangeWidth'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y<=self.filters[thisFiltName]
      
      # Chad's Matlab: % Filter out solutions with peaks wide in Doppler
      #   loop(:,26) - loop(:,25) <= maxDoppWidth
      # This looks like a type, as loop(:,25) is PKR2. He has a commented line:
      #   loop(:,27) - loop(:,26) <= maxDoppWidth
      # this would correspond to PKD2 - PKD1, which makes more sense for "Doppler Width".
      y = np.array(self.get_column('PKD2')) - np.array(self.get_column('PKD1'))
    
      thisFiltName = 'maxDoppWidth'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = y<=self.filters[thisFiltName]
      
      # Chad's Matlab: % Filter out solutions near Bragg
      # ~(abs(loop(:,9))>braggvel-maxcurrent & abs(loop(:,9))<braggvel+maxcurrent)
      # i.e., TRGV (target velocity) again. He uses
      # % Some Simple Calculations
      #  braggwavelength = 300/(2*txfreq);
      #  braggvel = sqrt(9.8*braggwavelength/(2*pi));  % Bragg Velocity (m/s)
      braggWavelength = 300./(2*self.TransmitCenterFreqMHz)
      braggVel = np.sqrt(9.8*braggWavelength/(2*np.pi)) # Bragg Velocity (m/s)
      y = np.abs(np.array(self.get_column('TRGV'))) # Target speed [m/s]
      thisFiltName = 'nearBragg'
      doKeep[thisFiltName] = np.array([True]*len(y))
      if ~np.isnan(self.filters[thisFiltName]):
        doKeep[thisFiltName] = np.logical_or(y<=(braggVel-self.filters[thisFiltName]),y>=(braggVel+self.filters[thisFiltName]))
      
      # doKeep for bearing limits--start with all true.
      y = np.array(self.get_column('TRGB')) # Target bearing
      doKeep['targetBearing'] = np.array([True]*len(y))
      for i in range(np.size(self.filters['bearingStart'])):
        bearingStart = self.filters['bearingStart'][i]
        bearingEnd = self.filters['bearingEnd'][i]
        isInSector = self.insector(bearingStart,bearingEnd,y)
        doKeep['targetBearing'] = np.logical_and(doKeep['targetBearing'],isInSector)
    
      # Combine all the indices to create a master index.
      #for filt in ('minSNR','maxSNR','minRange','maxRange','minVesselSpeed','maxVesselSpeed','minRangeWidth','maxRangeWidth','maxDoppWidth','nearBragg','targetBearing'):
      #  print("filt is",filt, "shape is",np.shape(doKeep[filt]))

      doKeep['master'] = np.array([True]*len(y))
      for filt in ('minSNR','maxSNR','minRange','maxRange','minVesselSpeed','maxVesselSpeed','minRangeWidth','maxRangeWidth','maxDoppWidth','nearBragg','targetBearing'):
        #print("filt is",filt)
        #print(np.shape(doKeep['master']),doKeep['master'])
        #print(np.shape(doKeep[filt]),doKeep[filt])
        doKeep['master'] = np.logical_and(doKeep['master'],doKeep[filt])
      
      print("filter_loop complete...")
      self.doKeep = doKeep

    ################################################################
    def insector(self,startBearing,endBearing,bearing,isInclusive=[True,True]):
        # insector--Returns numpy Boolean array which is True when numpy array of bearings (degrees)
        # lies between startBearing and endBearing. bearings can also be a scalar.
        #
        # Optional input argument isInclusive is [True,True] by default, i.e., endpoints are included. Change it to
        # [True,False] to exclude the endBearing, etc.
        #
        # Based on solution at https://stackoverflow.com/questions/22658954/check-if-angle-is-between-angles-from-and-to-with-clockwise-direction
        #
        # isInSector = insector(self,startBearing,endBearing,bearing,isInclusive)                                                                                                                                                       
        # Have to tell python not to go changing my variables of its own accord. 
        endBearing = np.copy(endBearing)
        bearing = np.copy(bearing)
        
        if not type(bearing) is np.ndarray:
            bearing = np.array(bearing)
        
        # Ensure that all inputs are between 0 and 360 degrees.
        if startBearing<0 or startBearing>=360:
            raise CodarLoopError("startBearing must be greater than or equal to 0 and less than 360 degrees.")
    
        if endBearing<0 or endBearing>=360:
            raise CodarLoopError("endBearing must be greater than or equal to 0 and less than 360 degrees.")
    
        if np.any(bearing<0) or np.any(bearing>=360):
            raise CodarLoopError("All elements of bearing must be greater than or equal to 0 and less than 360 degrees.")
    
        # If the end bearing is less than the start bearing, wrap it around 360 degrees.
        if endBearing <= startBearing:
            endBearing = endBearing + 360
    
        # Wrap bearing by 360 degrees wherever it is less than startBearing.
        bearing[bearing<startBearing] = bearing[bearing<startBearing] + 360
    
        if isInclusive[0]:
          isAbove = bearing>=startBearing
        else:
          isAbove = bearing>startBearing
    
        if isInclusive[1]:
          isBelow = bearing<=endBearing
        else:
          isBelow = bearing<endBearing
            
        isInSector = np.logical_and (isAbove,isBelow)
        return isInSector

if __name__ == '__main__':
    print("Running testdrive example")
    inFile = "./sampleData/vion/keep/LOOP_VION_2017_10_25_0000.loop.trunc"
    myLoop = Codar_loop(inFile)
    loopTime = myLoop.time[0]
    print("first loopTime is",loopTime[0])
    print("File date is",loopTime[0].date())
    print("TransmitCenterFreqMHz is ",myLoop.TransmitCenterFreqMHz)
    print("testdrive complete")

