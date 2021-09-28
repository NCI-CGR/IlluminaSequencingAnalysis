SS=`grep Error /DCEG/Projects/Exome/SequencingData/secondary_buf/logs/parse_build_familial_build_2019_23806/*.stdout | cut -d"." -f1 | cut -d"/" -f9`
for i in $SS; do
   SAMPLE=`echo $i | cut -d"_" -f2-`
   NUM=`grep NS5 /ttemp/wmy/${SAMPLE}_list | wc -l`
   echo
   if [[ $NUM -eq 4 ]]; then
     echo "found!"
   else
      cat /ttemp/wmy/${SAMPLE}_list
   fi 
done
