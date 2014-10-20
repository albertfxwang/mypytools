#!/usr/bin/env python

import os,glob

os.chdir('/home/khuang/mycode/madau')
rlist = glob.glob('tau*rest')
transdic = {}
wave = []
zs = []
print rlist[:3]

for i in range(len(rlist)):
    name = rlist[i]
    trans = []
    g = open(name)
    lines = g.readlines()
    for line in lines:
        data = line.split()
        if i == 0:
            wave += [float(data[0])]
        trans += [float(data[-1])]
    if i == 0:  # deal with duplicate wavelengths
        for j in range(len(wave)-1):  # wave is in decreasing order
            if wave[j] == wave[j+1]:
                wave[j] += 0.0001
        wave = wave[::-1]
    # figure our z
    z = name[4:7]
    z = float(z)
    zint = int(z*10)
    zs += [zint]
    transdic[zint] = trans[::-1]
# write to data file
f = open('madau_transmission.dat','w')
f.write('# 1 lam\n')
# write headers
for i in range(len(zs)):
    zk = zs[i]
    if zk < 10:
        f.write('# %d z0%d\n' % (i+2,zk))
    else:
        f.write('# %d z%d\n' % (i+2,zk))
# write data
for i in range(len(wave)):
    f.write('%f ' % wave[i])
    for j in range(len(zs)):
        zk = zs[j]
        f.write('%f ' % transdic[zk][i])
    f.write('\n')
f.flush()
f.close()
                
