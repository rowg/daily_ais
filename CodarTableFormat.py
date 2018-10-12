# Run testdrive with 
#   python3 CodarTableFormat.py

'''
REVISION HISTORY

2018-09-12, kpb--Backwards compatibility for Python 2 file i/o.

2018-10-05, kpb--Added test for truncated data file.

'''

class CodarTableFormatError(Exception):
   """Basic exception for CodarTableFormat class"""

class CodarTableFormat(object):
   # Class for reading CODAR Table Format (CTF) files. Table numbers sensibly 
   # start at 1 in original files, but are zero-indexed here because python.
   #
   # No reason to call as a script, so no input arguments to module. Will be 
   # called only from python; __init__() takes filename as argument.
   #
   # kpb@uvic.ca, 2017-10-31.
   def __init__(self,inFile,userDefinedColNames=None):
     # Call to constructor requires the name of the file to be read.
     self.inFile = inFile
     # Optional input argument is a dictionary object containing the column
     # names from some or all of the tables in the file. Specifying column
     # names over-rides the normal behaviour of using the column names in the
     # TableColumnTypes lines in the data file. userDefinedColNames is a
     # dictionary object with the keys being the zero-indexed table numbers.
     # For example,
     #   userDefinedColNames[0] = ['A', 'B', 'C']
     #   userDefinedColNames[2] = ['D', 'E', 'F', 'G']
     # (note that table #1 is not specified in this example--column names for
     # table 1 will be taken from table #1's TableColumnTypes line).
     self.userDefinedColNames = userDefinedColNames
     self.read_ctf()

   def calc_time(self, tableNum=None):
     print("calc_time() should be overloaded with a method written for the specific CTF file type (.ruv, .loop, etc.), since time format in these files can vary.")
     if tableNum is None:
       tableNum = self.find_primary_table()
     ctfTime = []
     return ctfTime
   
   def find_primary_table(self):
     # Returns the number of the primary table.
     for thisTableNum in self.data['tables']:
       if self.data['tables'][thisTableNum]['meta']['isPrimaryTable']:
         primaryTableNum = thisTableNum
         primaryTableFound = True
         break
     if not primaryTableFound:
       raise CodarTableFormatError("No primary table was found.")
     return primaryTableNum
     
   def get_column(self, colName, tableNum=None):
     # Get the named column. If tableNum is not specified, then the primary 
     # table is assumed.
     if tableNum is None:
       tableNum = self.find_primary_table()
       
     if tableNum not in self.data['tables']:
       raise CodarTableFormatError("Table specified by table number " + tableNum + " does not exist.")

     if colName not in self.data['tables'][tableNum]['data']:
       raise CodarTableFormatError("Data column specified by column name " + colName + " does not exist.")
       
     return self.data['tables'][tableNum]['data'][colName]

   def set_column(self, colName, x, tableNum=None):
     # Set the named column to the value x.

     # If tableNum is not specified, then the primary table is assumed.
     if tableNum is None:
       tableNum = self.find_primary_table()
       
     if tableNum not in self.data['tables']:
       raise CodarTableFormatError("Table specified by table number " + tableNum + " does not exist.")

     if colName not in self.data['tables'][tableNum]['data']:
       raise CodarTableFormatError("Data column specified by column name " + colName + " does not exist.")
       
     self.data['tables'][tableNum]['data'][colName] = x
     
   def read_ctf(self):
     # Reads in a CODAR Table Format (CTF) file into a python variable.
     # Warning: the reference file "File_CodarTableFormat.pdf" from CODAR Ocean
     # Sensors says that the TableColumnTypes (column names) entry for each 
     # table is optional--if omitted, then "the data's column order is assumed by
     # the type of file. This routine, however, requires its presence. May 
     # modify code later to permit the calling routine to specify column order... 
     #print("Call to read_ctf(), filename is ",self.inFile)
     self.data = {}
     self.data['fileMetaData'] = {}
     self.data['fileMetaData']['filename'] = self.inFile
     self.data['fileHeader'] = {}
     self.data['fileHeader']['variables'] = {}
     self.data['fileHeader']['comments'] = []
     self.data['fileFooter'] = {}
     self.data['fileFooter']['variables'] = {}
     self.data['fileFooter']['comments'] = []
     self.data['tables'] = {}
     self.data['time'] = {}

     try:
       inFid = open(self.inFile,'r',encoding = "ISO-8859-1")
     except TypeError:
       import io
       inFid = io.open(self.inFile,encoding = "ISO-8859-1")
     except Exception as e:
       print("Unexpected error. Failed to open loop file.")
       raise e

     tableNum = -1
     isInFileHeader = True
     isInTable = False
     isInTableHeader = False
     isInFileFooter = False
     fileEndFound = False

     line = 'dummy'
     while len(line)>0:
       line = inFid.readline().strip()
       if line.strip() == "%End:":
          fileEndFound = True

       # Detect whether this is a transition line into or out of a table or table header.
       if line[0:11] == "%TableType:":
         # Transition line into a table header.
         isInTableHeader = True
         tableNum = tableNum + 1

         # Initialise the new table.
         self.data['tables'][tableNum] = {}
         self.data['tables'][tableNum]['meta'] = {}
         self.data['tables'][tableNum]['tableHeader'] = {}
         self.data['tables'][tableNum]['tableHeader']['variables'] = {}
         self.data['tables'][tableNum]['tableHeader']['comments'] = []
         self.data['tables'][tableNum]['data'] = {}
         self.data['tables'][tableNum]['comments'] = []

         # ...Assume "primary", non-secondary table. Will be changed to False if table entries preceded by '% '.
         self.data['tables'][tableNum]['meta']['isPrimaryTable'] = True
         isInFileHeader = False

       elif line[0:12] == "%TableStart:":
         # Transition line into a table.
         isInTable = True
         isInTableHeader = False

         # Reading in of table header is complete, so can finish initialisation of table now by adding column names.
         # Use the TableColumnTypes CTF variable unless the column names have been specified in the constructor.
         doUseTableColumnTypes = True
         if self.userDefinedColNames is not None:
           if tableNum in self.userDefinedColNames.keys():
             # A list of column names has been specified for this table; use these names
             # instead of TableColumnTypes.
             doUseTableColumnTypes = False

         if doUseTableColumnTypes:
           colNames = self.data['tables'][tableNum]['tableHeader']['variables']['TableColumnTypes'].split()
         else:
           colNames = self.userDefinedColNames[tableNum]
             
         for colName in colNames:
           self.data['tables'][tableNum]['data'][colName] = []

       elif line[0:10] == "%TableEnd:":
         # Transition line leaving a table.
         isInTable = False

       else:
         # Not a transition line into/out of table, table header, etc.
         if isInFileHeader:
             if line[0:2] == "%%":
                 # This is a file header comment line.
                 comment = line[2:] # Leading "%%" stripped off comment
                 self.data['fileHeader']['comments'].append(comment)
             else:
                 # This is a file header data line.
                 cols = line.split(':')
                 varName = cols[0][1:] # Leading "%" stripped off variable name
                 varVal = cols[1]
                 self.data['fileHeader']['variables'][varName] = varVal

         elif isInTableHeader:
           if line[0:2] == "%%":
             # This is a table header comment line.
             comment = line[2:] # Leading "%%" stripped off comment
             self.data['tableHeader']['comments'].append(comment)
           else:
             # This is a table header data line.
             cols = line.split(':')
             varName = cols[0][1:] # Leading "%" stripped off variable name
             varVal = cols[1]
             self.data['tables'][tableNum]['tableHeader']['variables'][varName] = varVal

         elif isInTable:
           if line[0:2] == "%%":
             # This is a table comment line.
             # Table comment lines are too free-form to parse reliably. For example, .ruv files
             # have column headers with three rows, not all of which have the same number of
             # columns. Some, like "Radial Type", have no third line because there aren't
             # associated units; some have spaces in their names, etc.
             comment = line[2:] # Leading "%%" stripped off comment
             self.data['tables'][tableNum]['comments'].append(comment)
           else:
             # This is a table data line.
             if line[0:2] == "% ":
               # A table data line that starts with a '% ' denotes a table that is "secondary".
               self.data['tables'][tableNum]['meta']['isPrimaryTable'] = False
               # Strip off the leading '%' from the data line.
               line = line[1:]
             # The data in this line must now be read into the table-column variables.
             for iCol,colVal in enumerate(line.split()):
               thisColName = colNames[iCol]
               # Convert numeric values from string representation.
               try:
                 colVal = float(colVal)
               except:
                 # Actual strings (cannot be converted to float) are saved as-is.
                 pass
               self.data['tables'][tableNum]['data'][thisColName].append(colVal)

         else:
           # If not a transition line and not in file header, table header, or
           # table, then this is a file footer line.
           if line[0:2] == "%%":
             comment = line[2:] # Leading "%%" stripped off comment
             self.data['fileFooter']['comments'].append(comment)          
           elif line[0:5] != "%End:":
             cols = line.split(':')
             if len(cols) > 1:
               # Line is not an empty "%%" line.
               varName = cols[0][1:] # Leading "%" stripped off variable name
               varVal = cols[1]
               self.data['fileFooter']['variables'][varName] = varVal

     inFid.close()

     if not fileEndFound:
       raise CodarTableFormatError("No file end statement found. File is truncated/corrupted.")

if __name__ == '__main__':
    # Run testdrive by reading in a .ruv file. In practice, .ruv files will 
    # be read in using an ruv-specific subclass of the CodarTableFormat class.
    inFile = "./sampleData/vion/RDLm_VION_2017_10_23_0000.ruv"
    print("Running testdrive examples with inFile = " + inFile)

    # Create CodarTableFormat object using the default data column names 
    # (taken from the TableColumnTypes encoded in the file)
    CTF = CodarTableFormat(inFile)
    lond = CTF.get_column('LOND')
    print("First longitude (default variable name 'LOND') is",lond[0])

    # Create CodarTableFormat object using user-defined column names, rather
    # than the default ones in the file's TableColumnTypes.
    userDefinedColNames = {}
    userDefinedColNames[0] = ['longitude', 'latitude','U', 'V','vectorFlag', 'spatialQuality','temporalQuality', 'velocityMax','velocityMin', 'spatialCount','temporalCount', 'xDist','yDist', 'range', 'bearing','velocity','direction','rangeCell']
    CTF_userDefined = CodarTableFormat(inFile,userDefinedColNames)
    longitude = CTF_userDefined.get_column('longitude')
    print("First longitude (more descriptive, user-defined variable name 'longitude') is",longitude[0])

    # Get some data that is not in a data column (i.e., accessed directly,
    # rather than via get_column()). This works with .ruv and .loop files. Others?
    TransmitCenterFreqMHz = CTF_userDefined.data['fileHeader']['variables']['TransmitCenterFreqMHz']
    print("fileHeader.variables.TransmitCenterFreqMHz is ",TransmitCenterFreqMHz," MHz")

    # Print out the names of the available fields of the data dictionary.
    print("Fields of 'data' are: ", list(CTF_userDefined.data.keys()))
    print( "testdrive complete")
    

