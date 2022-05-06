#!/bin/env python
from __future__ import print_function
import os,sys

if __name__ == "__main__":
    search_term = sys.argv[1] 

    fecal_path = os.environ["FHICL_FILE_PATH"]
    #fecal_path = os.environ["FW_SEARCH_PATH"]
    fecal_dirs = fecal_path.split(":")
    num_current_search = 0
    for d in fecal_dirs:
        if d in [".","./"]:
            if num_current_search>0:
                continue
            else:
                num_current_search += 1
        if os.path.exists(d)==False:
            #print("Did not find directory: ",d)
            continue
        files = os.listdir(d)
        fecal_files = []
        for f in files:
            if ".fcl" in f:
                f1 = "%s/%s" % ( d, f )
                fecal_files.append(  f1 )
                #print(f1)
            if ".fcl" in search_term and f==search_term:
                print(f1)
            if os.path.basename(f.strip())==search_term:
                print("Found file: ",f.strip()," in ",d)
        for f in fecal_files:
            cmd = "grep -iIr %s %s" % ( search_term, f )
            pgrep = os.popen( cmd )
            lgrep = pgrep.readlines()
            if len(lgrep)>0:
                outname = 0
                for l in lgrep:
                    if "#include" in l:
                        continue
                    if outname==0:
                        print(cmd)
                        outname += 1
                    print(l)
