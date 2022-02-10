1. File description
   - pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.2021.csv: original Keytable
   - USU_CGR_ID_Table.csv: USU ID & CGR ID keytable
2. This delivery contains:
   - 638 normal subjets 
   - 6 topoffs 

3. You you can find the raw data in GCP:
   - https://console.cloud.google.com/storage/browser/dceg-covnet-wgs-useast1;tab=objects?project=nih-nci-dceg-cgr&prefix=&forceOnObjectsSortingFiltering=false

4. Original Keytable:
   - /home/lix33/lxwg/project/COVID19/02_10_2022/pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.2021.csv

5. Commandline used to generate this SampleID_Mapping_Table
```
$ tail -n +2 ./pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.2021.csv | awk -F ',' '{print $3}' | awk -F '-' '{print $NF}' > 3rdCol.txt
$ tail -n +2 ./pI3.COVNET.batch3.N638.6topoffs.keytable.dec.16.2021.csv | awk -F ',' '{print $1","$2}' | less > 1st2ndCol.txt
$ paste 1st2ndCol.txt 3rdCol.txt -d ',' > USU_CGR_ID_Table.csv
```

6. Data screenshot
![image](https://user-images.githubusercontent.com/11053933/153458093-35e802dd-28b7-4f51-8083-050cc2d7b8d4.png)
