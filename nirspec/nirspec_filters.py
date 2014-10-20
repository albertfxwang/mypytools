
def nirspec_filters(filter):
        disp = None
        smooth = None
	if filter=='NIRSPEC-1':
		start = 9470.
		end = 11210.
#                disp = (end-start)/800.
                disp = 2.15
                smooth = 20.
	elif filter=='NIRSPEC-2':
		start = 10890.
		end = 12970.
                disp = 2.13
                smooth = 18.
	elif filter=='NIRSPEC-3':
#		start = 11430.
#		end = 13750.
                start = 11000.
                end = 14250.
                disp = 2.79
	elif filter=='NIRSPEC-4':
		start = 12410.
		end = 15930.
	elif filter=='NIRSPEC-5':
		start = 14130.
		end = 18080.
		disp = 3.56
	elif filter=='NIRSPEC-6':
		start = 15580.
		end = 23150.
	elif filter=='NIRSPEC-7':
		start = 18390.
		end = 26300.
                disp = 4.08
        if disp is None:
           	disp = (start+end)/(2264.*4.)
        if smooth is None:
                smooth = 25.

	return start,end,disp,smooth
