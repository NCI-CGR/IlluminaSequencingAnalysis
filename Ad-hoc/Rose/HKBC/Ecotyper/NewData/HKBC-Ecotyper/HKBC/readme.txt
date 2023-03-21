(1) How to get tumor only columns
 
awk '
    BEGIN { FS=OFS="\t" }
    NR==1 {
        for (inFldNr=1; inFldNr<=NF; inFldNr++) {
            if ($inFldNr ~ /_T/) {
                out2inFldNr[++numOutFlds] = inFldNr
            }
        }
    }
    {
        for (outFldNr=1; outFldNr<=numOutFlds; outFldNr++) {
            inFldNr = out2inFldNr[outFldNr]
            printf "%s%s", $inFldNr, (outFldNr<numOutFlds ? OFS : ORS)
        }
    }' ./all.231tumor.normal.rsem.genes.results.tpm.txt | head -n 1 | tr '\t' '\n' | wc -l
