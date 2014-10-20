#!/usr/bin/env python

import os, glob, sys

"""
Rename the TFIT output files to include the specific run IDs.
"""

def rename_tfitrun(run_id):
   files = glob.glob('*')
   for f in files:
      root, ext = os.path.splitext(f)
      if not root.endswith(run_id):
         newname = root + '_' + run_id + ext
         os.system('mv %s %s' % (f, newname))

if __name__ == "__main__":
   run_id = sys.argv[1]
   rename_tfitrun(run_id)