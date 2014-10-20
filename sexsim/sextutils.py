# Utilities for parseing SExtractor files
#
# H. Ferguson - revised 4/23/03 to promote ints to floats if a value
# with a decimal point appears somewhere in the column originally thought
# to be integers
# 
# v2.1 - fails gracefully when the catalog has no sources
# v3.0 - added gettypes to return column types
#      - create new column names when they are not explictly in the header
# v4.0 - added gettypes to return column types
# v4.1 - uses numarray by default
# v4.2 - delete attributed 'l' (input lines from catalog) before returning
# v4.3 - 1/11/06Added less-offensive alias se_catalog() == sextractor()
# v4.4h- 1/21/06 Fixed bug in creating extra column names when last is a vector
# v4.4 - V. Laidler added new methods:
#           __len__ returns number of objects in catalog
#           __iter__ returns the index of the next row in the catalog
#           line(self,i) returns a constructed string containing the ith line 
#           buildheader returns a constructed header from the hdict
#        Added new attribute self.header: contains the header as read in
#           from the catalog.
#        Lines that start with '#' but are not followed by an integer are
#           now assumed to be comment lines, which are added to the
#           header but otherwise skipped.
# v4.5 - V. Laidler removed Numeric dependence
# v4.6 - V. Laidler converted to numpy
# v5.0 - 7/5/07 Numpy conversion
# v6.0 - V. Laidler: added rw_catalog class, reworked internals to avoid
#           column name clashes
# v7.0 - H. Ferguson: rewrote to allow appending rows and columns
#      - removed rw_catalog class, since se_catalog class now allows writing
#      - Redid the column type checking to also figure out the print format

__version__ = '6.0'
__author = 'Henry C. Ferguson, STScI'

import string
import numpy as N
import os, sys
import re

class se_catalog(object):
    """ Read a SExtractor-style catalog. 
        Usage: c=se_catalog(catalog,readfile=True,preserve_case=False,checkallformats=True)
        Will read the catalog and return an object c, whose attributes are 
        arrays containing the data. For example, c.mag_auto contains the 
        mag_auto values.
        Arguments:
        catalog -- The input SExtractor catalog. 
        readfile -- True means read the data. False means return the 
           object without reading the data. The lines from the catalog
           are returned as a list of ascii strings c.l. Useful if you want
           to do some special parsing of some sort. 
        preserve_case -- default (False) converts column names to lower case

        The input catalog MUST have a header with the SExtractor format:
           # 1 ID comment
           # 2 ALPHA_J200 another comment
        That is, first column is the comment symbol #, second column is
        the column number, third column is the column name, and the rest
        of the line is a comment. SExtractor allows "vectors" to be identified
        only by the first column...e.g.
           # 12 FLUX_APER
           # 20 FLUXERR_APER 
        the missing columns are all aperture fluxes through different
        apertures. These will be read into attributes:
           c.flux_aper   # The first one
           c.flux_aper_1 # the second one, and so on

        The case of aperture radii is a bit nasty, since these only
        appear in the SExtractor configuration file. Use parseconfig()
        to read that file.
    """
    def __init__(self,cfile,readfile=True,preserve_case=False,ncheck=100):
        # Initialize the catalog
        (d,lines,ncol,header) = initcat(cfile, preserve_case=preserve_case)
        # Save these as hidden attributes (so as not to confuse with columns)
        self._d = d # This is a dictionary of just the column names
        self._l = lines # All of the data rows of the catalog, as a list
        self._ncolumns = ncol # The number of columns
        self._header = header # The header
        self._fname = cfile # The catalog filename
        self._colnames=[] # The column names as an ordered list (for iteration)
        coldict = invert_dict(self._d)
        for k in coldict:
            self._colnames.append(coldict[k])
        if readfile:
            nlines = len(lines) # Number of data rows
            # Turn each line into a list, for faster access later
            self._colentries = range(nlines) 
            for i in range(nlines):
                self._colentries[i] = lines[i].split()
            # Check the formatting of at least some of the rows
            if ncheck == None or ncheck == 'all':  
                 ncheck = nlines
            self.gettypes(ncheck) # Only check a subset of the, for speed
            # Extract the columns
            for k in self._d.keys():
                contents = getcolvalues(self._d[k],self._type[k],
                                        self._colentries)
                try:
                    #Munge column name if it conflicts
                    test=self.__getattribute__(k)
                    newkey='c_'+k
                    print "--Column '%s' read in as '%s' to avoid conflicts"%(k,newkey)
                    setattr(self,newkey,contents)
                    self._d[newkey]=self._d[k]
                    del self._d[k]
                except AttributeError:
                    setattr(self,k,contents)
            delattr(self,'_l')

    def __len__(self):
        """ Return the number of rows """
        col = getattr(self,self._colnames[0]) # Grab a column 
        return len(col)                       # return its length

    def __iter__(self):
        return range(self.__len__).__iter__()

    def line(self,i):
        """ Returns an assembled line of this catalog suitable for writing.
        Except it doesn't really, if we modified the individual columns..."""
        ans = '    '.join(self._colentries[i])+'\n'
        return ans

    def buildheader(self):
        """ Reconstruct the header from the header dictionary.
        This might be useful if only a few columns were selected
        from the file; otherwise just use the 'header' attribute. """
        
        lines = {}
        for k in self._d:
            lines[self._d[k]]='#   %d   %s'%(self._d[k],k.upper())
        #sort the new keys
        nkeys= lines.keys()
        nkeys.sort()
        #join them together with newlines
        ans = ''
        for k in nkeys:
            ans=ans+"%s\n"%lines[k]
        return ans
        
    def gettypes(self,nrows):
        """ Figure out the type & format of each column from up to 100 rows """
        self._type = {}
        self._fmt = {}
        n  = len(self._l)
        step = max(1,n/nrows)
        for k in self._d.keys():
           ret = extract_column(self._d[k],self._l[:n:step]) # skip rows for speed
           fmt,tp,t = type_and_fmt(ret)
           self._type[k] = t
           self._fmt[k] = fmt
           
    def addcolumn(self, colname, coldata, fmt, comment=""):
        """ coldata must be a 1d numarray of the correct length"""
        if len(coldata) != len(self):
            raise ValueError,"Column length must match catalog length"

        #Most of the bookkeeping is the same as for an empty column
        self.addemptycolumn(colname,coldata.dtype,comment)

        #and then we reset the column to contain the actual data
        setattr(self,colname,coldata)
        self._fmt[colname] = fmt

    def addemptycolumn(self, colname, coltype, comment):
        """ Defines a new column & updates all the bookkeeping, but
        does not actually fill in the data. """
        setattr(self,colname,N.zeros((len(self),),coltype))
        self._type[colname] = 's'
        if coltype == 'float64':
            self._type[colname] = 'f'
        if coltype == 'int32':
            self._type[colname] = 'd'

        #Looks strange here because we count columns from 1 but
        #Python counts them from 0
        self._ncolumns += 1
        self._d[colname]=self._ncolumns
        self._colnames.append(colname)
        self._header += '# %d %s   %s\n'%(self._ncolumns,colname,comment)

    def row(self,rownum):
        """ Construct a new row to be printed out """
        line = ''
        for c in self._colnames:
            col=getattr(self,c)
#           print c, rownum, self._fmt[c], col[rownum] 
            line += self._fmt[c] % (col[rownum]) + ' '
        line += '\n'
        return line

    def writeto(self,outname,clobber=False):
        if not clobber:
            if os.path.isfile(outname):
                raise ValueError, """File already exists.
                   Use .writeto(fname, clobber=True) to overwrite. """

        out=open(outname,'w')

        out.write(self._header)
        for k in range(len(self)):
            out.write(self.row(k))
        out.close()

    def printme(self):
        """ Like writeto, but for sys.stdout """
        sys.stdout.write(self._header)
        for k in range(len(self)):
            sys.stdout.write(self.row(k))

class sextractor(se_catalog): # Just an alias for class se_catalog
    """ Read SExtractor catalog...just an alias for se_catalog """
    pass
class rw_catalog(se_catalog): # Just an alias for class se_catalog
    """ Read SExtractor catalog...just an alias for se_catalog """
    pass
def invert_dict(d):
    """ Generate a new dictionary with the key/value relationship inverted """
    newd={}
    for k in d:
        newd[d[k]]=k
    return newd

def parseconfig_se(cfile):
  """ parseconfig -- read a SExtractor .sex file and return a dictionary
    of options & values. Comments are ignored.
  """
  cdict = {}
  f = open(cfile,'r')
  lines = f.readlines()
  for l in lines:
    a = string.split(l)
    if len(a) > 0:
      if a[0][0] != '#':
        maxi = len(a)
        for i in range(1,len(a)):
          if a[i][0] == '#':
            maxi = i
            break
        # Turn comma-separated lists into python lists
        entry = []
        for e in a[1:maxi]:
          if string.find(e,','):
            entry = entry + string.split(e,',')
          else:
            entry = entry + [e]
        cdict[a[0]] = entry
  return cdict

def initcat(catfile,preserve_case=False):
  """ parseheader -- reads the header of a SExtractor catalog file and 
      returns a dictionary of parameter names and column numbers.
      Also returns a list of lines containing the data.
  """
  hdict = {}
  header=[]
  f = open(catfile,'r')
  lines = f.readlines()
  f.close()
  first = 1
  firstdata = 0
  i = 0
  previous_column=0
  previous_key=""
  for l in lines:
    if l.startswith('#'): #this is a header line        
        header.append(l)        
        a = (l.replace('#','# ')).split() #Guard against "#10 colname"
        try:
            col = int(a[1])
            # If the column numbers skip, create new column names for
            # columns not named explicitly in the header
            if col != previous_column+1:
              for c in range(previous_column+1,col):
                column_name = previous_key+"_%d" % (c-previous_column)
                hdict[column_name] = c
            # Update this column in the dictionary
            if (preserve_case):
                column_name = a[2]
            else:
                column_name = a[2].lower()
                hdict[column_name] = col
            firstdata = i+1
            previous_column = col
            previous_key = column_name
        except (ValueError, IndexError):
            #it's a comment line with no column number,
            #or an entirely blank comment line: skip
            pass
    
    else:  # This is where the data start
        a=string.split(l)
        if len(a)>0:
            if first:
              firstdata = i
              first = 0
            # Check if there are extra columns
            if len(a) > previous_column: 
            # If so, add keys for the last entry
              for c in range(previous_column+1,len(a)+1):
                column_name = previous_key+"_%d" % (c-previous_column)
                if (preserve_case):
                    hdict[column_name] = c
                else:
                    hdict[column_name.lower()] = c
            ncolumns = len(a)
    i = i+1
  return(hdict,lines[firstdata:],ncolumns,''.join(header))

def getcolvalues(col,coltype,colentries):
  """ Get a column from a SExtractor catalog. Determine the type
      (integer, float, string) and return either an array of that
      type (Int32, Float64) or a list of strings """
  i = col-1               # Columns start at 1, arrays start at 0
  nlines = len(colentries)
  if len(colentries) == 0:
    values = N.array([])
    return values
  a = colentries	
  
  if coltype == 's':	  # Convert strings
      values = range(nlines)
      for j in range(nlines):
          values[j] = colentries[j][i]
      values = N.array(values)

  if coltype == 'f':	# Convert floats 
      values = N.zeros((nlines),N.float64)
      for j in range(nlines):
          values[j] = float(colentries[j][i])

  if coltype == 'd':	# Convert Ints
      values = N.zeros((nlines),N.int32)
      for j in range(nlines):
          values[j] = int(colentries[j][i])
  return values

def writeheader(fh,colnames):
    """ Write an SExtractor-style header to an open file handle.

    @param fh: file handle
    @type fh: file

    @param colnames: list of column names
    @type colnames: list

    @todo: add space checking to colnames
    @todo: permit passing a filename?
    @todo: handle comments
    """
    for i in range(len(colnames)):
        fh.write('#   %d   %s\n'%(i+1,colnames[i]))

def type_and_fmt(stringlist):
    """ Determine the type and format of a list of strings.
        The longest string determines the length and the most digits after
	the decimal point determines the precision. 
	Precedence: 
	   - strings
	   - exponential format floats
	   - floats
	   - ints
    """
    var_type = type('')
    length = 0
    floatlen = 0
    stringlen = 0 
    precision = 0
    fmt_type = 'd'
    for l in stringlist:
        f = fmt(l)
        stringlen = max(stringlen,f.l)
	if f.fmt_type == fmt_type:
            length = max(length,f.l)
            precision = max(precision,f.p)
	if f.fmt_type == 's':
            fmt_type = 's'
            length = max(length,f.l)
	    precision = 0
	    var_type = f.t
            length = max(length,f.l)
	if f.fmt_type == 'e' and fmt_type != 's':
            # Reset the length and precision if this is the first instance
            if fmt_type != 'e':
                savelen = max(length,f.l)
	        length = f.l 
	        precision = f.p
            fmt_type = 'e'
	    var_type = f.t
	    length = f.l
	    precision = max(precision,f.p)
	if f.fmt_type == 'f' and fmt_type != 's' and fmt_type != 'e':
          val = float(l)
          if not N.isinf(val) and not N.isnan(val):
            if fmt_type != 'f': # Looks like we are changing from ints to floats on this row!
                floatlen = max(length,len(str(int(float(l))))) # max number of digits before the decimal pt.
	        precision = f.p # number of digits after the decimal point
                length = max(floatlen+precision+1,f.l)
            else: 
                floatlen = max(floatlen,len(str(int(float(l))))) # max number of digits before the decimal pt.
	        precision = max(precision,f.p) # number of digits after the decimal point
                length = max(floatlen+precision+1,f.l)
            fmt_type = 'f'
	    var_type = f.t
    if fmt_type == 's':
        outfmt = "%%%ds" % (stringlen)
    if fmt_type == 'd':
        outfmt = "%%%dd" % (length)
    if fmt_type == 'f':
        outfmt = "%%%d.%d%s" % (length,precision,fmt_type)
    # Once we've used the 'e' for the format, we can tell the caller it's 'f'
    if fmt_type == 'e':  
        outfmt = "%%%d.%d%s" % (length+2,precision,fmt_type) # To allow space for e+00 exponent
        fmt_type = 'f'
    return (outfmt, var_type,fmt_type)
        
class fmt:
    """ Determine the format of a string 
        fmt.fmt = the format sfor printing the string
	fmt.t = the type for converting to float, string, or int
	fmt.l = the length of the string
	fmt.precision = the precision of floats
    """
    def __init__(self,s):	 
        # Output format is going to be %<length>.<precision><fmt_type>
        length = len(s)
        precision = 0
        fmt_type = ''
    
        # Check if this can be converted to an integer
        try:
           int(s)
           stype = type(1)
           fmt_type = 'd'	 
        except:
            # If not, check if it can be converted to a float
            try: 
                float(s)
            except:
                fmt_type = 's' # If not, assume it's a string
    	        stype = type('')
	if fmt_type != 'd' and fmt_type != 's':
            fmt_type = 'f'
    	    stype = type(1.2)
            ss = re.split('[eE]',s) # See if this is in exponentiated format
            if len(ss) > 1:
                fmt_type = 'e'
            sss = ss[0].split('.')  # Find the decimal point
            if len(sss) > 1:
                precision = len(sss[1]) # number of digits after decimal pt.
        # Here is the derived format
	if precision == 0:
            self.fmt = "%%%d%s" % (length,fmt_type)
	else:
            self.fmt = "%%%d.%d%s" % (length,precision,fmt_type)
        self.t = stype
	self.l = length
	self.p = precision
	self.fmt_type = fmt_type

def extract_column(col,lines):
    """ Extract a column from a bunch of rows """
    ret = []
    i = col-1
    for l in lines:
       a = string.split(l)[i]
       ret = ret + [a]
    return ret
