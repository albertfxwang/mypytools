#!/bin/tcsh -m 
#=======================================================================
#align HST data to Subaru
#=======================================================================
#+
# NAME:
#   align_HSTSubaru.csh 
#
# PURPOSE:
#   correct coordinates in HST catalog to match Subaru images
#
# COMMENTS:
#
#
# INPUTS:
#   HST image Subaru image catalog              two catalogues
#
#   --Subaru   Subaru image
#   --HST Hst image
#   --cat catalog with entries ALPHA_J2000 and DELTA_J2000  
# OPTIONAL INPUTS:
#   --nosex No running of sextractor: You better be sure you have 
#	    all the ref catalogs though 
# OUTPUTS:
#   new catalog
#   new HST image with corrected WCS
#
# EXAMPLES:
#   
# BUGS:
#   #  
# REVISION HISTORY:
#-
#=======================================================================

set simage = 0
set himage = 0
set cat = 0
set help = 0
set sex = 1
set sexall = 1
set match = 0.001
set sexconf=0
set telescope = 0
while ( $#argv > 0 )
    switch ($argv[1])
	case --{help}:   
	      shift argv
	      set help = 1
	      breaksw
        case --{Subaru}:
	     shift argv  	
	     set simage = $argv[1]
             shift argv  
           breaksw
        case --{HST}:
	     shift argv  
	     set himage = $argv[1]
	     set telescope = 1
             shift argv  
           breaksw
        case --{Spitzer}:
	     shift argv  
	     set himage = $argv[1]
	     set telescope = 2
             shift argv  
           breaksw
       case --{nosex}:
	     shift argv  
	     set sex = 0 
           breaksw
       case --{match}:                
             shift argv                                                    
             set match = $argv[1]
             shift argv
           breaksw                                                            
       case --{matchedcat}:                
             shift argv                                                    
             set merged = $argv[1]
	     set sexall = 0
	     set sex = 0
             shift argv
           breaksw    
        case --{sexconf}:                
             shift argv                                                    
             set sexconf = $argv[1]
             shift argv
           breaksw    
                                                       
         case *: 
            set cluster = $argv[1]
            shift argv
            breaksw
        endsw
end
if ($help) then
  print_script_header.csh $0
  goto FINISH
endif


 if ((! -e $simage) || (! -e $himage)) then
        echo "We need Subaru image $simage or HST image  $himage " ; exit 1
 endif
if ($sexall == 1) then 
 
 set refcat = HSTSUB.cat
if ($sexconf == 0) then
 set sexconf = $SHAGGLES_DIR/sex/shag_detect_noweight.sex
endif
source $SHAGGLES_DIR/bin/progs_csh.ini
 set dir = ./
 if ($sex == 1) then
    if ((! -e $sexconf) || (! -e $himage)) then
	    echo "Need $sexconf configuration file" ; exit 1
    endif
    echo "Running sextractor"


    ${P_SEX} -c $sexconf -CATALOG_NAME ${refcat:r}.hscat  ${himage} -MAG_ZEROPOINT 25.80 

    ${P_SEX} -c $sexconf -SATUR_LEVEL 50000.0 -CATALOG_NAME ${refcat:r}.sscat  ${simage} -MAG_ZEROPOINT 25.80 

    echo "switching to $dir"
    cd $dir 
    
    sex2imcat_marusa.pl -r -o ${refcat:r}.hcat ${refcat:r}.hscat 
    sex2imcat_marusa.pl -r -o ${refcat:r}.scat ${refcat:r}.sscat 
 
    lc 'x = %ALPHA_J2000 %DELTA_J2000 2 vector'  < ${refcat:r}.scat > tmp1.cat
    lc 'x = %ALPHA_J2000 %DELTA_J2000 2 vector'  <  ${refcat:r}.hcat > tmp2.cat
    mkreg.pl -wcs -ds9 -ds9clean -xcol 0 -ycol 1 -colour red -rad 2 tmp1.cat
    mkreg.pl -wcs -ds9 -xcol 0 -ycol 1 -colour green -rad 1 tmp2.cat
    echo "Red is subaru, green is hst"
 endif #sex =1
 mergecats $match tmp1.cat tmp2.cat | lc > merged.imcat
 set nm = `lc -c < merged.imcat`
 echo "$nm objects when merged"
 mkreg.pl -ds9 -wcs -xcol 0 -ycol 1 -colour blue -rad 3 merged.imcat
 echo "Blue is original that were found"
 lc -o < merged.imcat > merged.cat   
endif 
set xmin = `catstats -s < merged.imcat | lc -o 'x = %x[0][0]' | head -1`
set xmax = `catstats -s < merged.imcat | lc -o 'x = %x[0][0]' | head -2 | tail -1`
set ymin = `catstats -s < merged.imcat | lc -o 'x = %x[0][1]' | head -1`
set ymax = `catstats -s < merged.imcat | lc -o 'x = %x[0][1]' | head -2 | tail -1`
 
echo "$xmin $xmax $ymin $ymax"
set himageout = ${himage:r}_aligned_py.fits
touch geo.tmp
echo "from myiraf import *" > geo.tmp
echo 'mygeotran("merged.cat", "'"$himage"'", "'"$himageout"'"'",$xmin,$xmax,$ymin,$ymax, fitgeometry='shift')" >> geo.tmp
python < geo.tmp

 set xshift = `more geomap.db | grep xshift | cut -f 4`
 set yshift = `more geomap.db | grep yshift | cut -f 4`
 #set rot = `more geomap.db | grep xrotation | cut -f 3`
 #set rotn=`awk 'BEGIN{print (('$rot' > 180) ? '-360+$rot' : '$rot')}'`
set rotn = 0
set tolmax = '0.005'
set tolmin = '-0.005'

set converged=`awk 'BEGIN{print (('$xshift' > '$tolmax') || ('$xshift' < '$tolmin') || ('$yshift' >  '$tolmax')  || ('$yshift' < '$tolmin')) ? 0 : 1}'`

if ($converged == 0) then
 echo "Shifts very big $xshift $yshift $rotn"
endif
 
set pixscale = `imsize $himage | cut -d '/' -f 2 | cut -d 's' -f 1`


set xshiftp=`echo $xshift $pixscale | awk '{printf "%10.12f\n", -$1*(3600.0*cos($1*3.14/180.0))/$2;}'`
set yshiftp=`echo $yshift $pixscale | awk '{printf "%10.12f\n", $1*3600.0/$2;}'`


echo "Shifts $xshiftp $yshiftp $rotn (pixscale $pixscale)"

if ($telescope == 1) then
    ur_setup
    set himageout = ${himage:r}_aligned.fits
    /bin/cp ${himage} ${himageout}
    set  hdrinput = "align.hd.input"
    set hdroutput = "align.hd.output"
    echo "Updating align.hd.input"
    if ( -f $hdrinput ) \rm -f $hdrinput
    if ( -f $hdroutput ) \rm -f $hdroutput
    echo "from drizzlepac import updatehdr" >! $hdrinput
    echo "from stwcs import updatewcs" >> $hdrinput
    #echo "updatewcs.updatewcs('*flt.fits')" >> $hdrinput
    echo "updatehdr.updatewcs_with_shift('${himageout}', '${himage}',rot=${rotn}, scale=1.0, xsh=${xshiftp}, ysh=${yshiftp}, verbose=True, force=True)" >> $hdrinput
    python ${hdrinput}
else if ($telescope == 2) then
    set himageout = ${himage:r}_aligned.fits
    /bin/cp ${himage} ${himageout}
     echo "Updating $himageout"
    set crval1 = `gethead $himageout CRVAL1`
    set crval2 = `gethead $himageout CRVAL2`
    set crval1n=`echo $xshift $crval1 | awk '{printf "%10.12f\n", $1 + $2;}'`
    set crval2n=`echo $yshift $crval2 | awk '{printf "%10.12f\n", $1 + $2;}'`
    echo $crval1 $crval2
    echo $crval1n $crval2n
    keyhead -r $himageout CRVAL1=$crval1n
    keyhead -r $himageout CRVAL2=$crval2n
else
    echo "Don't know which telescope I'm dealing with!"
endif

FINISH:

#=======================================================================
