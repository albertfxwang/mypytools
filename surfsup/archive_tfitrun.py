#!/usr/bin/env python

import os, sys
import rename_tfitrun

def archive(run_id):
   if not os.path.exists(run_id):
      os.mkdir(run_id)
   os.system('mv *collage*.fits %s/' % run_id)
   os.system('mv *.cat_best* %s/' % run_id)
   os.system('mv ddiags* %s/' % run_id)
   os.system('mv tpipe*.log %s/' % run_id)
   os.system('mv *.cell %s/' % run_id)
   os.system('mv dlog* %s/' % run_id)
   #os.system('cp *.param %s/' % run_id)
   #os.system('cp ../sub_medbkgd_JD.py %s/' % run_id)
   os.chdir(run_id)
   rename_tfitrun.rename_tfitrun(run_id)

if __name__ == "__main__":
   run_id = sys.argv[1]
   archive(run_id)
   