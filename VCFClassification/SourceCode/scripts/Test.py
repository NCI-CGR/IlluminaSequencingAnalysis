import os
outVCF = "/home/lixin/lxwg/Data/ad-hoc/VCFClusification/OSTEOSARCOMA.genesOfInterest.rare.nonsyn.vcf"
outVCFTmp = outVCF + ".tmp"
outVCFOrg = outVCF + ".org"
if os.path.exists(outVCF):
    if os.path.exists(outVCFTmp):
        CMD = "rm " + outVCFTmp
        os.system(CMD)
    if os.path.exists(outVCFOrg):
        CMD = "rm " + outVCFOrg
        os.system(CMD)
    CMD = "sed '/#CHROM/q' " + outVCF + " >> " + outVCFTmp
    os.system(CMD)
    CMD = "sed -n '/#CHROM/,$p' " + outVCF + " | tail -n +2 " + "| sort -k 2,2 | uniq >> " + outVCFTmp
    os.system(CMD)
    CMD = "mv " + outVCF + " " + outVCFOrg
    os.system(CMD)
    CMD = "mv " + outVCFTmp + " " + outVCF
    os.system(CMD)
    
    
     