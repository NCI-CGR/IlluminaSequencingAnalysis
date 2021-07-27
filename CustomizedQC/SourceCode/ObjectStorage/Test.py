'''
Move data from biowulf to s3 
This is a simple testing, including dry_run and wet_run
'''
import subprocess
import sys
import os
from traitlets.config.application import catch_config_error

def main():
    strDataDir = "/home/lix33/lxwg/Test/object_storage"
    strRootDir = "/home/lix33/lxwg/Test"
    strObjPrefix = "lix33/Test/ObjStorage/DirTest"
    
    #Get all of files 
    CMD = "find " + strDataDir + " -type f"
    strFastqList = subprocess.getoutput(CMD)
    vFile = strFastqList.split('\n')
    for strFile in vFile:
        if strFile == "":
            continue
        strPrefix = strObjPrefix + os.path.dirname(strFile).replace(strRootDir, '') + "/"        
        CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " --dry-run -V" # Dry run
        #CMD = "obj_put -v DCEG_COVID_WGS -p " + "\"" + strPrefix + "\" " + strFile + " -V" # wet run
        print(CMD)
        iReturn = os.system(CMD)
        if iReturn != 0:
            print("\n", ">>>>>> Error Occured!")
            print(CMD, "\n", "<<<<<<<", "\n") 
            exit
    print("All Set!")

if __name__ == "__main__":    
    main()
    