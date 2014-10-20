"""
Calculates the helio-centric velocity of a given source at a given time. This
  is just a python transcription of some other code (I forgot who's code I
  used...).
"""

from math import pi,cos,sin,atan,asin,sqrt

AUKM = 1.4959787e8
CTROP = 365.24219572
CBES = 0.313
C1900 = 1900.0
C2000 = 2415020.0
DEGRAD = pi/180.

def bcv_keck(mjd,ra,dec,eq):
	return bcv(mjd,155.478333,19.8283333333333,4160.,ra,dec,eq)

def bcv(mjd,long,lat,alt,ra,dec,eq):
	radlong = long*DEGRAD
	radlat = lat*DEGRAD

	radra = ra*DEGRAD
	raddec = dec*DEGRAD

	st = sidereal_time(mjd,radlong)

	eqt = (mjd-C2000-CBES)/CTROP + C1900
	c = [cos(radra)*cos(raddec),sin(radra)*cos(raddec),sin(raddec)]
	prema = pre(eq,eqt)

	cc = []
	for i in range(3):
		cc.append(c[0]*prema[i][0]+c[1]*prema[i][1]+c[2]*prema[i][2])

	if cc[0]!=0:
		arg = cc[1]/cc[0]
		ra2 = atan(arg)
		if cc[0]<0:
			ra2 += pi
		elif cc[1]<0:
			ra2 += 2.*pi
	else:
		if cc[1]>0:
			ra2 = pi/2.
		else:
			ra2 = 1.5*pi

	dec2 = asin(cc[2])

	ha = st-ra2

	gcvel = geovel(radlat,alt,dec2,-1.*ha)
	velh,velb = barvel(mjd,eqt)

	bcvel = 0.
	hcvel = 0.

	for i in range(3):
		bcvel += velb[i]*cc[i]*AUKM
		hcvel += velh[i]*cc[i]*AUKM

	return gcvel,bcvel,hcvel

def sidereal_time(mjd,deglong):
	DTPI = 2.0*pi
	DF = 1.00273790934
	DCT0 = 2415020.0
	DCJUL = 36525.0

	D1 = 1.739935934667999
	D2 = 6.283319509909095e2
	D3 = 6.755878646261384e-6

	jd0 = int(mjd) + 0.5
	if jd0>mjd:
		jd0 -= 1.
	ut = (mjd-jd0)*DTPI
	t = (jd0-DCT0)/DCJUL
	st0 = D1 + D2*t + D3*t*t
	st0 = st0%DTPI
	st = DF*ut + st0 - deglong
	st = (st+2*DTPI)%DTPI

	return st

def pre(eq1,eq2):
	DCSAR = 4.848136812e-6
	DC1900 = 1900.0
	DC1M2 = 0.01
	DC1 = 2304.25
	DC2 = 1.396
	DC3 = 0.302
	DC4 = 0.018
	DC5 = 0.791
	DC6 = 2004.683
	DC7 = -0.853
	DC8 = -0.426
	DC9 = -0.042

	DT0 = (eq1 - DC1900)*DC1M2
	DT = (eq2 - eq1)*DC1M2
	DTS = DT * DT
	DTC = DTS * DT
	DZETA = ((DC1+DC2*DT0)*DT+DC3*DTS+DC4*DTC)*DCSAR
	DZETT = DZETA + DC5*DTS*DCSAR
	DTHET = ((DC6+DC7*DT0)*DT+DC8*DTS+DC9*DTC)*DCSAR
	DSZETA = sin(DZETA)
	DCZETA = cos(DZETA)
	DSZETT = sin(DZETT)
	DCZETT = cos(DZETT)
	DSTHET = sin(DTHET)
	DCTHET = cos(DTHET)
	DA = DSZETA * DSZETT
	DB = DCZETA * DSZETT
	DC = DSZETA * DCZETT
	DD = DCZETA * DCZETT

	dprema = []
	dprema.append([DD*DCTHET-DA,-1.*DC*DCTHET-DB,-1.*DSTHET*DCZETT])
	dprema.append([DB*DCTHET+DC,-1.*DA*DCTHET+DD,-1.*DSTHET*DSZETT])
	dprema.append([DCZETA*DSTHET,-1.*DSZETA*DSTHET,DCTHET])

	return dprema

def geovel(DPHI,DH,DEC,DHA):
	DA = 6378.140
	DF = 0.00335281
	DW = 7.2921158554e-5

	DE2 = DF*(2.-DF)

	D1 = 1.-DE2*(2.-DE2)*sin(DPHI)**2
	D2 = 1.-DE2*sin(DPHI)**2
	DR0 = DA*sqrt(D1/D2)

	D1 = DE2*sin(2.*DPHI)
	D2 = 2.*D2
	DPHIG = DPHI-atan(D1/D2)

	DRH = DR0*cos(DPHIG) + (DH/1000.)*cos(DPHI)
	DVELG = DW*DRH*cos(DEC)*sin(DHA)

	return DVELG

def barvel(DJE,DEQ):

	DC2PI = 6.2831853071796
	CC2PI = 6.283185
	DC1 = 1.
	DCT0 = 2415020.0
	DCJUL = 36525.0
	DCBES = 0.313
	DCTROP = 365.24219572
	DC1900 = 1900.0

	DCFEL = [[1.7400353,6.2833195099091e2,5.2796e-6],
			[6.2565836, 6.2830194572674e2,-2.6180e-6],
			[4.7199666, 8.3997091449254e3,-1.9780e-5],
			[1.9636505e-1, 8.4334662911720e3,-5.6044e-5],
			[4.1547339, 5.2993466764997e1, 5.8845e-6],
			[4.6524223, 2.1354275911213e1, 5.6797e-6],
			[4.2620486, 7.5025342197656, 5.5317e-6],
			[1.4740694, 3.8377331909193, 5.6093e-6]]

	DCEPS = [4.093198e-1,-2.271110e-4,-2.860401e-8]

	CCSEL = [[1.675104e-2,-4.179579e-5,-1.260516e-7],
			[2.220221e-1, 2.809917e-2, 1.852532e-5],
			[1.589963, 3.418075e-2, 1.430200e-5],
			[2.994089, 2.590824e-2, 4.155840e-6],
			[8.155457e-1, 2.486352e-2, 6.836840e-6],
			[1.735614, 1.763719e-2, 6.370440e-6],
			[1.968564, 1.524020e-2,-2.517152e-6],
			[1.282417, 8.703393e-3, 2.289292e-5],
			[2.280820, 1.918010e-2, 4.484520e-6],
			[4.833473e-2, 1.641773e-4,-4.654200e-7],
			[5.589232e-2,-3.455092e-4,-7.388560e-7],
			[4.634443e-2,-2.658234e-5, 7.757000e-8],
			[8.997041e-3, 6.329728e-6,-1.939256e-9],
			[2.284178e-2,-9.941590e-5, 6.787400e-8],
			[4.350267e-2,-6.839749e-5,-2.714956e-7],
			[1.348204e-2, 1.091504e-5, 6.903760e-7],
			[3.106570e-2,-1.665665e-4,-1.590188e-7]]

	DCARGS = [[5.0974222,-7.8604195454652e2],
			 [3.9584962,-5.7533848094674e2],
			 [1.6338070,-1.1506769618935e3],
			 [2.5487111,-3.9302097727326e2],
			 [4.9255514,-5.8849265665348e2],
			 [1.3363463,-5.5076098609303e2],
			 [1.6072053,-5.2237501616674e2],
			 [1.3629480,-1.1790629318198e3],
			 [5.5657014,-1.0977134971135e3],
			 [5.0708205,-1.5774000881978e2],
			 [3.9318944, 5.2963464780000e1],
			 [4.8989497, 3.9809289073258e1],
			 [1.3097446, 7.7540959633708e1],
			 [3.5147141, 7.9618578146517e1],
			 [3.5413158,-5.4868336758022e2]]

	CCAMPS = [[-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5,-2.490817e-7],
			[-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5,-1.823138e-7],
			[6.593466e-7, 1.322572e-5, 9.258695e-6,-4.674248e-7,-3.646275e-7],
			[1.140767e-5,-2.049792e-5,-4.747930e-6,-2.638763e-6,-1.245408e-7],
			[9.516893e-6,-2.748894e-6,-1.319381e-6,-4.549908e-6,-1.864821e-7],
			[7.310990e-6,-1.924710e-6,-8.772849e-7,-3.334143e-6,-1.745256e-7],
			[-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6,-1.655307e-7],
			[-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6,-3.736225e-7],
			[3.442177e-7, 2.671323e-6, 1.832858e-6,-2.394688e-7,-3.478444e-7],
			[8.702406e-6,-8.421214e-6,-1.372341e-6,-1.455234e-6,-4.998479e-8],
			[-1.488378e-6,-1.251789e-5, 5.226868e-7,-2.049301e-7, 0.],
			[-8.043059e-6,-2.991300e-6, 1.473654e-7,-3.154542e-7, 0.],
			[3.699128e-6,-3.316126e-6, 2.901257e-7, 3.407826e-7, 0.],
			[2.550120e-6,-1.241123e-6, 9.901116e-8, 2.210482e-7, 0.],
			[-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.]]

	CCSEC3 = -7.757020e-8

	CCSEC = [[1.289600e-6, 5.550147e-1, 2.076942],
			[3.102810e-5, 4.035027, 3.525565e-1],
			[9.124190e-6, 9.990265e-1, 2.622706],
			[9.793240e-7, 5.508259, 1.559103e1]]

	DCSLD = 1.990987e-7
	CCSGD = 1.990969e-7

	CCKM = 3.122140e-5
	CCMLD = 2.661699e-6
	CCFDI = 2.399485e-7

	DCARGM = [[5.1679830, 8.3286911095275e3],
			[5.4913150,-7.2140632838100e3],
			[5.9598530, 1.5542754389685e4]]

	CCAMPM = [[1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7],
			[-2.223581e-2, 5.083103e-8, 1.002548e-2,-2.291823e-8],
			[1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8]]

	CCPAMV = [8.326827e-11,1.843484e-11,1.988712e-12,1.881276e-12]

	DC1MME = 0.99999696

	FORBEL = []
	SORBEL = []
	SN = [0.,0.,0.,0.]
	SINLP = [0.,0.,0.,0.]
	COSLP = [0.,0.,0.,0.]
	for i in range(7):
		FORBEL.append(float(i))
	for i in range(17):
		SORBEL.append(float(i))

	IDEQ = DEQ
	DT = (DJE - DCT0)/DCJUL
	T = DT
	DTSQ = DT * DT
	TSQ = DTSQ

	for k in range(8):
		DLOCAL = (DCFEL[k][0]+DT*DCFEL[k][1]+DTSQ*DCFEL[k][2])%DC2PI
		if k==0:
			DML = DLOCAL
		else:
			FORBEL[k-1] = DLOCAL
	DEPS = (DCEPS[0]+DT*DCEPS[1]+DTSQ*DCEPS[2])%DC2PI

	for k in range(17):
		SORBEL[k] = (CCSEL[k][0]+T*CCSEL[k][1]+TSQ*CCSEL[k][2])%DC2PI

	for k in range(4):
		A = (CCSEC[k][1]+T*CCSEC[k][2])%DC2PI
		SN[k] = sin(A)

	PERTL = CCSEC[0][0]*SN[0] + CCSEC[1][0]*SN[1] + (CCSEC[2][0]+T*CCSEC3)*SN[2] + CCSEC[3][0]*SN[3]
	PERTLD = 0.
	PERTR = 0.
	PERTRD = 0.

	for k in range(15):
		A = (DCARGS[k][0]+DT*DCARGS[k][1])%DC2PI
		COSA = cos(A)
		SINA = sin(A)
		PERTL += CCAMPS[k][0]*COSA + CCAMPS[k][1]*SINA
		PERTR += CCAMPS[k][2]*COSA + CCAMPS[k][3]*SINA
		if k>=10:
			continue
		PERTLD += (CCAMPS[k][1]*COSA-CCAMPS[k][0]*SINA)*CCAMPS[k][4]
		PERTRD += (CCAMPS[k][3]*COSA-CCAMPS[k][2]*SINA)*CCAMPS[k][4]

	E = SORBEL[0]
	G = FORBEL[0]

	ESQ = E**2
	DPARAM = DC1 - ESQ
	PARAM = DPARAM

	TWOE = E*2.
	TWOG = G*2.

	PHI = TWOE*((1.-ESQ*(1./8.))*sin(G) + E*(5./8.)*sin(TWOG) + ESQ*0.5416667*sin(G+TWOG))
	F = G+PHI
	SINF = sin(F)
	COSF = cos(F)
	DPSI = DPARAM/(DC1+E*COSF)
	PHID = TWOE*CCSGD*((1.+ESQ*1.5)*COSF+E*(1.25-SINF*SINF*0.5))
	PSID = CCSGD*E*SINF/sqrt(PARAM)

	D1PDRO = (DC1 + PERTR)
	DRD = D1PDRO*(PSID+DPSI*PERTRD)
	DRLD = D1PDRO*DPSI*(DCSLD+PHID+PERTLD)
	DTL = (DML+PHI+PERTL)%DC2PI
	DSINLS = sin(DTL)
	DCOSLS = cos(DTL)
	DXHD = DRD*DCOSLS - DRLD*DSINLS
	DYHD = DRD*DSINLS + DRLD*DCOSLS

	PERTL = 0.
	PERTLD = 0.
	PERTP = 0.
	PERTPD = 0.

	for k in range(3):
		A = (DCARGM[k][0] + DT*DCARGM[k][1])%DC2PI
		SINA = sin(A)
		COSA = cos(A)
		PERTL   = PERTL  + CCAMPM[k][0]*SINA
		PERTLD  = PERTLD + CCAMPM[k][1]*COSA
		PERTP   = PERTP  + CCAMPM[k][2]*COSA
		PERTPD  = PERTPD - CCAMPM[k][3]*SINA

	TL =  FORBEL[1] + PERTL
	SINLM = sin(TL)
	COSLM = cos(TL)
	SIGMA = CCKM/(1.0 + PERTP)
	A = SIGMA*(CCMLD+PERTLD)
	B = SIGMA*PERTPD
	DXHD = DXHD + A*SINLM + B*COSLM
	DYHD = DYHD - A*COSLM + B*SINLM
	DZHD = 0.0  - SIGMA*CCFDI*cos(FORBEL[2])

	DXBD = DXHD*DC1MME
	DYBD = DYHD*DC1MME
	DZBD = DZHD*DC1MME

	for k in range(4):
		PLON = FORBEL[k+3]
		POMG = SORBEL[k+1]
		PECC = SORBEL[k+9]

		TL = (PLON+2.*PECC*sin(PLON-POMG))%CC2PI
		SINLP[k] = sin(TL)
		COSLP[k] = cos(TL)

		DXBD += CCPAMV[k]*(SINLP[k]+PECC*sin(POMG))
		DYBD += -1.*CCPAMV[k]*(COSLP[k]+PECC*cos(POMG))
		DZBD += -1.*CCPAMV[k]*SORBEL[k+13]*cos(PLON-SORBEL[k+5])

	DCOSEP = cos(DEPS)
	DSINEP = sin(DEPS)
	DYAHD = DCOSEP*DYHD - DSINEP*DZHD
	DZAHD = DSINEP*DYHD + DCOSEP*DZHD
	DYABD = DCOSEP*DYBD - DSINEP*DZBD
	DZABD = DSINEP*DYBD + DCOSEP*DZBD

	if IDEQ==0:
		DVELH = [DXHD,DYAHD,DZAHD]
		DVELB = [DXBD,DYABD,DZABD]
		return DVELH,DVELB

	DEQDAT = (DJE - DCT0 - DCBES)/DCTROP + DC1900
	DPREMA = pre(DEQDAT,DEQ)

	DVELH = []
	DVELB = []
	for i in range(3):
		DVELH.append(DXHD*DPREMA[i][0]+DYAHD*DPREMA[i][1]+DZAHD*DPREMA[i][2])
		DVELB.append(DXBD*DPREMA[i][0]+DYABD*DPREMA[i][1]+DZABD*DPREMA[i][2])
	return DVELH,DVELB
