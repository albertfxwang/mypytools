#!/usr/bin/env python
import pydoc
import sys
import os.path
import os
import glob
import copy
import shutil

home = '/goodssoft/tools/pygoods/'

modules = [
'angsep.py',
'coords.py',
'parseconfig.py',
'numprint.py',
'readcol.py',
'sextutils.py',
'match.py '
]

if len(glob.glob('doctmp')) == 0:
    os.mkdir('doctmp')

for i in range(len(modules)-1):
    infile = home + modules[i]
    outfile = 'doctmp/'+modules[i]
    f = open(infile)
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
        if lines[i].find('import *') >= 0:
           lines[i] = '# ' + lines[i]
    f = open(outfile,'w')
    f.writelines(lines)

os.chdir('doctmp')
sys.path = [os.getcwd()] + sys.path
files = glob.glob('*.py')
for f in files:
    print f
    module = f[:-3]
    exec("import " + module)
    sys.stdout = open("../%s.txt" % (module),'w')
    s = help(module)
    sys.argv = ['pydoc','-w',module]
    pydoc.cli()
os.chdir('..')
shutil.rmtree('doctmp')
