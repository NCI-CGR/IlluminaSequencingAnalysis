package nci.cgr.manifest;

import java.io.IOException;
import java.text.ParseException;

public class SecondaryParser {

	public static void main(String[] args) throws IOException, InterruptedException, ParseException {
		// TODO Auto-generated method stub
//		for(int i = 0; i < args.length; i++) {
//            System.out.println(args[i]);
//        }
		// The changes comparing to first version:
	    // 1. Somatic or germline is passed to step2 script;
		// 2. scan all the fastq files with new quality trimming;
		//
		if (args.length==0) usage(); 
		else{
			boolean found=false;
			if (args[0].equals("parse")){
				//parse $MANIFEST ${BUFFER_DIR}/input/all_merging_records.txt 
				// ${OUTPUT_DIR}/step2_merge_${DATE}.sh ${OUTPUT_DIR}/manifest_errs.txt 
				// ${OUTPUT_DIR}/new_update_analysisIDs.txt ${OUTPUT_DIR}/filenames_restore.txt 
				// ${MANIFEST_UPDATE} $DP $REDO ${BUFFER_DIR}/input/qt_all.txt
				
				// e.g.
				// parse T:\DCEG\Projects\Exome\builds\build_familial_build_2019_23806\Manifest\FAMILIAL-PLUS-CONTROLS-MANIFEST-7.16.2019_unix_noratio.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\all_merging_records.txt T:\DCEG\Home\wangm6\tmp.sh T:\DCEG\Home\wangm6\tmp_errs.txt T:\DCEG\Home\wangm6\tmp_new_update.txt T:\DCEG\Home\wangm6\tmp_restore.txt T:\DCEG\Home\wangm6\tmp_manifest.txt 60 update T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\qt_all.txt
				
				
				if (args.length==11){
					found=true;
					Manifest manifest=new Manifest(args[1],args[2],args[9]);
					manifest.parseManifest(args[3], args[4],args[5],args[6],args[7],args[8],args[10]);
				}
				
			}
			if (args[0].equals("scan")){
				//scan T:\\DCEG\\Projects\\Exome\\SequencingData\\variant_scripts\\logs\\GATK T:\DCEG\Home\wangm6\tmp2\all_merging_records.txt T:\\DCEG\\Projects\\Exome\\SequencingData\\variant_scripts\\logs\\all_merging_err.txt 
				if (args.length==4) {
					found=true;
					LogMerging logMerging=new LogMerging(args[1]);
					logMerging.scanLogs(args[2],args[3]);
				}
			}
			
			if (args[0].equals("compare")){
				//compare T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\bam_results_build_familial_build_2019_23806\manifest.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\bam_results_build_familial_build_2019_23806\all_out.txt T:\DCEG\Home\wangm6\tmp2\report.txt T:\DCEG\Home\wangm6\tmp2\error.txt
				if (args.length==5) {
					   found=true;
					   BamInfo bamInfo=new BamInfo();
					   String fileRecords=args[1];
					   String fileBamOutput=args[2];
					   String fileReport=args[3];
					   String fileErr=args[4];
					   bamInfo.parseBamInfo(fileRecords,fileBamOutput,fileReport,fileErr);
				}
			}
		    if (args[0].equals("qt")){  // Compare current manifest to the qt_all file to see if have any lane-level samples still in old trimming
		    	                        // Print out the fastq files which are needed restoration and sample list for retrimming
		    	
		    	// e.g.
		    	// qt T:\DCEG\Projects\Exome\builds\build_IBMFS-Familial_2018_21413\Manifest\IBMFS-7-24-2018.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\qt_all.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\fastq_restoration.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\samples_retrimming.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\qt_err.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\samples_new_trimmed.txt
		    	// qt T:\DCEG\Projects\Exome\builds\test_new_PLCO\Manifest\SR0390-TESTING-MANIFEST.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\qt_all.txt H:\tmp\fastq_restoration.txt H:\tmp\samples_retrimming.txt H:\tmp\qt_err.txt H:\tmp\new_trimmed.txt
		    	if (args.length==7){
		    		found=true;
		    		String fileErr=args[5];
		    		Manifest manifest=new Manifest(args[1]);
		    		manifest.verifyManifest(fileErr);
		    		manifest.compareNewQualityTrimmingList(args[2], args[3], args[4],args[6],fileErr);
		    	}
		    }
		    if (!found)  usage();
		}
	}
	
	  static void usage() {
		  
		  // Test configuration 1:
		  // qt H:\\projects\secondary_analysis\\AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv T:\\DCEG\\Projects\\Exome\\SequencingData\\sync_script_v2\\input\\qt_all_20171217.txt T:\\DCEG\\Home\\wangm6\\tmp2\\fastq_restoration.txt T:\\DCEG\\Home\\wangm6\\tmp2\\samples_retrimming.txt T:\\DCEG\\Home\\wangm6\\tmp2\\qt_err.txt  T:\\DCEG\\Home\\wangm6\\tmp2\\samples_new_trimmed.txt
		  //
		  // Test configuration 2:
		  // parse T:\DCEG\Projects\Exome\builds\build_familial_build_2019_23806\Manifest\FAMILIAL-PLUS-CONTROLS-MANIFEST-7.16.2019_unix_noratio.csv 
		  // T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\all_merging_records.txt 
		  // T:\DCEG\Home\wangm6\tmp2\step2_merge.sh 
		  // T:\DCEG\Home\wangm6\tmp2\step2_err.sh 
		  // T:\DCEG\Home\wangm6\tmp2\new_update_analysis.txt 
		  // T:\DCEG\Home\wangm6\tmp2\file_restoration.txt 
		  // T:\DCEG\Home\wangm6\tmp2\FMAILIAL-POP-CONTROL_update.txt 
		  // 60 update T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\qt_all.txt 
		  
		    System.err.println("(Version 2.1) Usage:");
	        System.err.println("1. java -jar secondaryParsing.jar parse <Manifest file> <Existing merging record file> <Step2 shell script file> <Error file> <Update or new analysis ID file>"
	        		+"<File listing for restoration> <Updated manifest file> <Threshold of coverage> redo|update <New quality trimming file list>");
	        System.err.println("Example: java -jar secondaryParsing.jar parse /DCEG/Projects/Exome/SequencingData/Manifest/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv "
	        		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/logs/all_merging_records.txt /DCEG/Projects/Exome/SequencingData/sync_script_v2/step2_merge.sh "
	        		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/manifest_errs.txt /DCEG/Projects/Exome/SequencingData/sync_script_v2/output/new_update_analysisIDs.txt "
	        		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/filenames_restore.txt /DCEG/Projects/Exome/SequencingData/Manifest/AUG-FMAILIAL-POP-CONTROL_update.txt 60 redo germline"
	        		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/input/qt_all.txt");
	        System.err.println();
	        System.err.println("2. java -jar secondaryParsing.jar scan <LOG path> <Merging record output file> <Error file>");
		    System.err.println("Example: java -jar secondaryParsing.jar scan /DCEG/Projects/Exome/SequencingData/variant_scripts/logs/GATK "
		    		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/input/all_merging_records.txt /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/input/all_merging_err.txt");
	        
		    System.err.println("3. java -jar secondaryParsing.jar compare <Merging record output file> <BAM record output file> <Report File> <Error file>");
		    System.err.println("Example: java -jar secondaryParsing.jar compare /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/input/all_merging_records.txt "
		    		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/bam_results/all_out.txt /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/cmp_report.txt"
		    		+" /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/cmp_err.txt");
	        // compare T:\DCEG\Projects\Exome\builds\build_SR0403-001_Mexican_breast_cancer_2019_22639\Manifest\NP0403-HE5_6_ANALYSIS_MANIFEST.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\bam_results_build_SR0403-001_Mexican_breast_cancer_2019_22639_2019-02-20\all_out.txt T:\DCEG\Home\wangm6\tmp2\cmp_report.txt T:\DCEG\Home\wangm6\tmp2\cmp_err.txt 
		    
		    System.err.println("3. java -jar secondaryParsing.jar qt <Manifest file> <Exisitng new quality trimming file list> <Fastq file restoration list> <Sample list for quality trimming redo>"
		    		+ " <Error file>");
		    System.err.println("Example: java -jar secondaryParsing.jar qt /DCEG/Projects/Exome/SequencingData/Manifest/AUG-FAMILIAL-POP-CONTROL-MANIFEST.csv /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/input/qt_all.txt "
		    		+"/DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/fastq_restoration.txt /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/samples_retrimming.txt"
		    		+" /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/qt_err.txt /DCEG/Projects/Exome/SequencingData/secondary_script_v2.1/output/samples_newtrimmed.txt");
		    // qt T:\DCEG\Projects\Exome\builds\build_IBMFS-Familial_2018_21413\Manifest\IBMFS-7-24-2018.csv T:\DCEG\Projects\Exome\SequencingData\secondary_buf\input\qt_all.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\fastq_restoration.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\samples_retrimming.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\qt_err.txt T:\DCEG\Projects\Exome\SequencingData\secondary_buf\output\20180809\samples_new_trimmed.txt
		    System.exit(-1);
	    }

}
