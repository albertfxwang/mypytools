# H. Ferguson 7/12/05
def degdms(ideg):
    if (ideg < 0):
       s = -1
    else:
       s = 1
    ideg = abs(ideg)
    deg = int(ideg)+0.
    m = 60.*(ideg-deg)
    minutes = int(m)+0.
    seconds = 60.*(m-minutes)
    if s < 0:
       dms = "-%02d:%02d:%06.3f" % (deg,minutes,seconds)
    else:
       dms = "%02d:%02d:%06.3f" % (deg,minutes,seconds)
    return dms

def deghms(ideg):
    ihours = ideg/15.
    hours = int(ihours)+0.
    m = 60.*(ihours-hours)
    minutes = int(m)+0.
    seconds = 60.*(m-minutes)
    hms = "%02d:%02d:%06.3f" % (hours,minutes,seconds)
    return hms

def to_sexagesimal(ra,dec):
    hms = deghms(ra)
    dms = degdms(dec)
    str = "%s %s" % (hms,dms)
    return(str)

def parsecoo(coostring,inhours=False):
    """ Convert an coordinate string to decimal degrees. Accept
        the input in any one of the following formats:
               hhh:mm:ss.s    -- where hhh is either hours or degrees
               hhh:mm.m       -- can include whitespace
               hhh mm ss.s
               hhh mm.m
               ddd.ddddd  -- assumed to be in degrees, not hours
               ###h##m##.###s -- with or without whitespace
               ###d##m##.###s -- with or without whitespace
        If inhours is True, multiply result by 15. (hours->deg)
    """
    if inhours:
        factor = 15.
    else: 
        factor = 1.0
    s = coostring.lstrip()  # Strip leading whitespace
    s = s.rstrip()         # Strip trailing whitespace
    s = s.replace(':',' ')     # Change colons to spaces
    s = s.replace('d',' ')     
    s = s.replace('D',' ')     
    s = s.replace('h',' ')     
    s = s.replace('H',' ')     
    s = s.replace('m',' ')     
    s = s.replace('M',' ')     
    s = s.replace('s',' ')     
    s = s.replace('S',' ')     
    s = s.replace('\'',' ')     
    s = s.replace('\"',' ')     
    a = s.split()          # Split the sting
#   print "coostring ", a
    if len(a) == 1:        # Assume degrees if just one value
        coo = float(a[0])
    if a[0][0] == '-':
        sgn = -1
        degorhours = -float(a[0])
    else:
        sgn = 1
        degorhours = float(a[0])
    if len(a) == 2:        # Assume hours, minutes if two values given
        coo = degorhours + float(a[1])/60.
        coo = sgn*factor*coo
    if len(a) == 3:        # Assume hours, minutes, seconds
        coo = degorhours + float(a[1])/60. + float(a[2])/3600.
        coo = sgn*factor*coo
    return coo

def parsera(rastring):
    return parsecoo(rastring,inhours=True)

def parsedec(decstring):
    return parsecoo(decstring)

def parsecoords(coostring):
    """ Convert coordinate string to decimal degrees. Accept
        the input in any one of the following formats:
        hhh:mm:ss.s ddd:mm:ss.s  -- whitespace allowed, as long as consistent
        hhh:mm:ss.s ddd:mm:ss.s  
        hhh mm ss.s ddd mm ss.s
        hhh:mm.m ddd:mm.m
        hhh mm.m ddd mm.m
        ddd.ddddd ddd.ddddd      -- ra,dec both assumed to be decimal deg
        
        also converts:
        12h25m32.52s 25d32m25.2s -- with or without whitespace
    """
    a = coostring.split()
    l = len(a)
    if l<2 or l%2 != 0:
        print "parsecoords: cannot convert %s" % coostring
        return
    if (len(a)) == 2:
        ra = parsera(a[0])
        dec = parsedec(a[1])
    if (len(a)) == 4:
        ra = parsera(a[0]+' '+a[1])
        dec = parsedec(a[2]+' '+a[3])
    if (len(a)) == 6:
        ra = parsera(a[0]+' '+a[1]+' '+a[2])
        dec = parsedec(a[3]+' '+a[4]+' '+a[5])
    return ra,dec
