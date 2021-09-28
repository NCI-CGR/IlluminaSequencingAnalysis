package nci.cgr.manifest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.ListIterator;
import java.nio.file.*;

import nci.cgr.manifest.*;

public class Manifest {

	public Manifest(String fileName,String fileRecords,String redo) {
		super();
		this.fileName = fileName;
		this.fileNameCurrentAllRecords=fileRecords;
		this.allExistingAnalysisIDs=null;
		this.newAnalysisIDs=null;
		this.redo=redo;
		// this.variantType=variantType.toUpperCase();
	}

	
	public Manifest(String fileName) {
		super();
		this.fileName = fileName;
	
	}
	
	public boolean verifyManifest(String fileErr) throws IOException{
		File file= new File(fileName);
		if (!file.exists()){
			System.err.println("Error: "+fileName+" is not existed!");
			return true;
		}
	    BufferedReader reader = new BufferedReader(new FileReader(file));
	    BufferedWriter writerErr=new BufferedWriter(new FileWriter(fileErr));
	    
	    String line="";
	    String analysisID="";
	    String flowcell="";
	    String lane="";
	    String sampleID="";
	    String index="";
	    String group="";
	    String instrumentID="";
	    String platform="";
	    String seqDate="";
	    Date seqDate2=null;
	    Boolean foundError=false;
	    float coverage=0;
	    System.out.println("Verifying the manifest file ...");
	    newAnalysisIDs = new ArrayList<LaneSampleList>();
	    line=reader.readLine();
	    String[] fields=line.split(",");
	    if (fields.length!=14){
	    	writerErr.append("Error: the header line does not have 14 columns!");
            writerErr.newLine();
	    }
	    int nodeCount=0;
	    int lineCount=0;
	    writerErr.append("# Parsing the manifest file: "+fileName);
	    writerErr.newLine();
	    while ((line=reader.readLine())!=null){
	    	String aa[]=line.split(",");
	    	lineCount++;
	    	
	    	// System.out.println(lineCount+":"+aa.length);
	    	if ((aa.length!=14) && (aa.length!=13)){
	 	       		int lineNum=lineCount+1;
	         		System.out.println("Error: Line number "+lineNum+" ("+line+") contains less than 14 or 13 columns!");
	    	    	writerErr.append("Error: Line number "+lineNum+" ("+line+") contains less than 14 or 13 columns!");
	    	        writerErr.newLine();
	    	        continue;
	    	    }
	     
	 	    
	    	// Manifest format requirement
	    	analysisID=aa[12];  //Analysis ID
	    	flowcell=aa[2];		// Flowcell 
	    	lane=aa[3];         // lane
	    	index=aa[4];
	    	sampleID=aa[5];
	    	group=aa[6];
	    	instrumentID=aa[0];
	    	seqDate=aa[1];
	    	if (analysisID.contains(" ") || flowcell.contains(" ") || lane.contains(" ") || index.contains(" ") || sampleID.contains(" ") || group.contains(" ") || instrumentID.contains(" ")){
	    		writerErr.append("Error: "+line+" contains space in fields!");
				writerErr.newLine();
	    	}
	    	
	    	 String expectedPattern = "MM/dd/yyyy";
	    	 SimpleDateFormat formatter = new SimpleDateFormat(expectedPattern);
	    	    
	    	 try
	    	 {
	    	      // (2) give the formatter a String that matches the SimpleDateFormat pattern
	    	     
	    	      seqDate2 = formatter.parse(seqDate);

	    	      // (3) prints out "Tue Sep 22 00:00:00 EDT 2009"
	    	      // System.out.println(date);
	    	  }
	    	  catch (ParseException e)
	    	    {
	    	      // execution will come here if the String that is given
	    	      // does not match the expected format.
	    	      e.printStackTrace();
	    	      writerErr.append("Error: "+line+" the date format is error!");
				  writerErr.newLine();
	    	    }
	    	    
	    	coverage=Float.parseFloat(aa[10]);
	    	boolean found=false;
	    	int i=0;
	    	for(i=0;i<newAnalysisIDs.size();i++){
	    	   if (newAnalysisIDs.get(i).analysisID.equals(analysisID)){
	    		   found=true;
	    		   break;
	    	   }
	    	}
	    	if (instrumentID.startsWith("E0131")|| instrumentID.startsWith("E0357")|| instrumentID.startsWith("E0047")){
	    		writerErr.append("Error: "+line+" the intrument ID is not HiSeq or NextSeq!");
				writerErr.newLine();
	    		foundError=true;
	    		continue;
	    	}
	    	else{
	    		if (instrumentID.startsWith("E0296"))
	    			platform="NextSeq";
	    		else if (instrumentID.startsWith("E0265"))
	    			platform="MiSeq";
	    		else
	    		    platform="HiSeq";
	    	}
	    	if (!found){
	    		LaneSampleNode sampleNode=new LaneSampleNode(flowcell, lane, sampleID,group, index,platform,"", "",0,0,"","","",1,coverage,null,seqDate2);
	    		boolean foundExistingError=false;
	    		for(i=0;i<newAnalysisIDs.size();i++){
	 	    	   LaneSampleNode existingNode= newAnalysisIDs.get(i).head;
	 	    	   while (existingNode!=null){
	 	    		  if (sampleNode.compareTo(existingNode)==true){
	 	    			 foundExistingError=true;
	 	    		     break;
	 	    		  }
	 	    		  existingNode=existingNode.next;
	 	    	   }
	 	    	   if (foundExistingError==true){
	 	    		    writerErr.append("Error: "+line+" is duplicate with other records!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	 	    	   }
	 	    		   
	 	    	}
	    		if (!foundError){
	    		   LaneSampleList analysisIDUnit=new LaneSampleList(analysisID,sampleNode,group);
	    		   newAnalysisIDs.add(analysisIDUnit);
	    		   nodeCount++;
	    		}
	    	}
	    	else{
	    		LaneSampleNode storedSampleNode=newAnalysisIDs.get(i).head;
	    		LaneSampleNode newSampleNode=new LaneSampleNode(flowcell, lane, sampleID,group, index, platform,"","",0,0,"","","",1,coverage,null,seqDate2);
	    		do
	    		{
	    			if (storedSampleNode.compareTo(newSampleNode)==true){
	    				writerErr.append("Error: "+line+" is duplicate with other records within the same analysis ID!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	    			}
	    			if (!storedSampleNode.group.equals(group)){
	    				writerErr.append("Error: "+line+" group name is different to other lane samples within one analysis ID!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	    			}
	    			if (!storedSampleNode.sampleID.equals(sampleID)){
	    				writerErr.append("Warning: "+line+" CGF ID is different to other lane samples within one analysis ID!");
	    				writerErr.newLine();
	    			}
	    			storedSampleNode=storedSampleNode.next;
	    		}
	    		while (storedSampleNode!=null);
	    		if (!foundError){
		    		newAnalysisIDs.get(i).add(newSampleNode);
		    		nodeCount++;
	    		}
	    	}
	    }
	    System.out.println("Verification is done! Total: "+nodeCount+" lines was finished processing!");
	    reader.close();
	    writerErr.append("# Manifest verification is done! Total: "+nodeCount+" lines");
		writerErr.newLine();
		writerErr.append("");
		writerErr.newLine();
	    writerErr.flush();
	    writerErr.close();
	   
		return foundError;
	}
	
	public void compareNewQualityTrimmingList(String fileNewQualityTrimmingList, String fileFastqFileRestorationList, String fileSampleListForNewQualityTrimming, String fileNewTrimmedSamples, String fileErr) throws IOException, ParseException, InterruptedException{
		System.out.println("Comparing to existing new quality trimming Fastq files ...");
		FastqFileList fastqFileList=new FastqFileList(fileNewQualityTrimmingList);
		BufferedWriter writerRestoration=new BufferedWriter(new FileWriter(fileFastqFileRestorationList));
		BufferedWriter writerSampleRedo=new BufferedWriter(new FileWriter(fileSampleListForNewQualityTrimming));
		BufferedWriter bwErr=new BufferedWriter(new FileWriter(fileErr,true));
		BufferedWriter writerSampleNewTrimmed=new BufferedWriter(new FileWriter(fileNewTrimmedSamples));  // The samples are already new trimmed
		int ret=fastqFileList.loadFastqRecords(bwErr);
		if (ret!=0)
			System.out.println("There are errors in FASTQ QT file list! Please check details in "+fileErr);
    	
    	for (int j=0;j<newAnalysisIDs.size();j++){
    		System.out.println("Comparing for "+j);
    		LaneSampleList sampleList=newAnalysisIDs.get(j);
  //  		if (sampleList.analysisID.equals("NA12878_GDNA_HS4K_KHPL_1"))
  //  			System.out.println("found!");
    		if (sampleList.compareToFastqFileList(fastqFileList,writerRestoration,writerSampleRedo)) {//all lanes are new trimmed
    			writerSampleNewTrimmed.append(sampleList.analysisID);
    			writerSampleNewTrimmed.newLine();
    			writerSampleNewTrimmed.flush();
    			System.out.println(sampleList.analysisID+" is new-trimmed!");
    		}	    		
     	}
    	writerRestoration.close();
		writerSampleRedo.close();
		writerSampleNewTrimmed.close();
		bwErr.close();
		System.out.println("Comparison is done!");
    	
    		
	}

	
	String fileName="";
	String fileNameCurrentAllRecords="";
	String redo="";
	// String dedup;
	// String variantType="GERMLINE";
	List<LaneSampleList> allExistingAnalysisIDs;
	List<LaneSampleList> newAnalysisIDs;
	
	public  void parseManifest(String fileOutput,String fileErr,String fileListUpdate,String fileNeedRestore,String fileManifestNew,String requestedDepth,String fileNewQualityTrimming) throws IOException, InterruptedException, ParseException {
		File file= new File(fileName);
		if (!file.exists()){
			System.err.println("Error: "+fileName+" is not existed!");
			return;
		}
		// System.out.println("Info: 2019-04-10");
		
	    BufferedReader reader = new BufferedReader(new FileReader(file));  //Read manifest csv file
	    BufferedWriter writerErr=new BufferedWriter(new FileWriter(fileErr));
	    BufferedWriter writerManifest=new BufferedWriter(new FileWriter(fileManifestNew));
	    
	    String line="";
	    String analysisID="";
	    String flowcell="";
	    String lane="";
	    String sampleID="";
	    String index="";
	    String group="";
	    String instrumentID="";
	    String platform="";
	    float coverage=0;
	    System.out.println("Verifying the manifest file and reformat as another manifest file ...");
	    newAnalysisIDs = new ArrayList<LaneSampleList>();
	    line=reader.readLine();
	    String[] fields=line.split(",");
	    writerManifest.append(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[5]+"\t"+fields[6]+"\t"+fields[7]+"\t" 
	    		+fields[12]+"\t"+fields[11]+"\t"+fields[13]+"\t"+fields[10]+"\t"+fields[8]+"\t"+fields[9]);
	    writerManifest.newLine();
	    
	    Boolean foundError=false;
	    int nodeCount=0;
	    int lineCount=0;
	    writerErr.append("# Parsing the manifest file: "+fileName);
	    writerErr.newLine();
	    while ((line=reader.readLine())!=null){
	    	String aa[]=line.split(",");
	    	lineCount++;
	    	// System.out.println(lineCount+":"+aa.length);
	    	if (aa.length==13){
	 	       writerManifest.append(aa[0]+"\t"+aa[1]+"\t"+aa[2]+"\t"+aa[3]+"\t"+aa[4]+"\t"+aa[5]+"\t"+aa[6]+"\t"+aa[7]+"\t" 
	 	    		+aa[12]+"\t"+aa[11]+"\t"+"\t"+aa[10]+"\t"+aa[8]+"\t"+aa[9]);
	 	        writerManifest.newLine();  
	    	}
	    	else{ 
	    		if (aa.length==14){
	    	       writerManifest.append(aa[0]+"\t"+aa[1]+"\t"+aa[2]+"\t"+aa[3]+"\t"+aa[4]+"\t"+aa[5]+"\t"+aa[6]+"\t"+aa[7]+"\t" 
	 	    		+aa[12]+"\t"+aa[11]+"\t"+aa[13]+"\t"+aa[10]+"\t"+aa[8]+"\t"+aa[9]);
	    	       writerManifest.newLine();
	    		}
	         	else{
	         		int lineNum=lineCount+1;
	         		System.out.println("Error: Line number "+lineNum+" ("+line+") contains less than 13 columns!");
	    	    	writerErr.append("Error: Line number "+lineNum+" ("+line+") contains less than 13 columns!");
	    	        writerErr.newLine();
	    	        continue;
	    	    }
	    	}
	 	    
	    	// Manifest format requirement
	    	analysisID=aa[12];  //Analysis ID
	    	flowcell=aa[2];		// Flowcell 
	    	lane=aa[3];         // lane
	    	index=aa[4];
	    	sampleID=aa[5];
	    	group=aa[6];
	    	instrumentID=aa[0];
	    	
	    	if (analysisID.contains(" ") || flowcell.contains(" ") || lane.contains(" ") || index.contains(" ") || sampleID.contains(" ") || group.contains(" ") || instrumentID.contains(" ")){
	    		writerErr.append("Error: "+line+" contains space in fields!");
				writerErr.newLine();
	    	}
	    	coverage=Float.parseFloat(aa[10]);
	    	boolean found=false;
	    	int i=0;
	    	for(i=0;i<newAnalysisIDs.size();i++){
	    	   if (newAnalysisIDs.get(i).analysisID.equals(analysisID)){
	    		   found=true;
	    		   break;
	    	   }
	    	}
	    	if (instrumentID.startsWith("E0131")|| instrumentID.startsWith("E0357")|| instrumentID.startsWith("E0047")){
	    		writerErr.append("Error: "+line+" the intrument ID is not HiSeq or NextSeq!");
				writerErr.newLine();
	    		foundError=true;
	    		continue;
	    	}
	    	else{
	    		if (instrumentID.startsWith("E0296"))
	    			platform="NextSeq";
	    		else if (instrumentID.startsWith("E0265"))
	    			platform="MiSeq";
	    		else
	    		    platform="HiSeq";
	    	}
	    	if (!found){
	    		LaneSampleNode sampleNode=new LaneSampleNode(flowcell, lane, sampleID,group, index,platform,"", "",0,0,"","","",1,coverage,null);
	    		boolean foundExistingError=false;
	    		for(i=0;i<newAnalysisIDs.size();i++){
	 	    	   LaneSampleNode existingNode= newAnalysisIDs.get(i).head;
	 	    	   while (existingNode!=null){
	 	    		  if (sampleNode.compareTo(existingNode)==true){
	 	    			 foundExistingError=true;
	 	    		     break;
	 	    		  }
	 	    		  existingNode=existingNode.next;
	 	    	   }
	 	    	   if (foundExistingError==true){
	 	    		    writerErr.append("Error: "+line+" is duplicate with other records (in flowcell, lane, index and sampleID) but in different analysis IDs!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	 	    	   }
	 	    		   
	 	    	}
	    		if (!foundError){
	    		   LaneSampleList analysisIDUnit=new LaneSampleList(analysisID,sampleNode,group);
	    		   newAnalysisIDs.add(analysisIDUnit);
	    		   nodeCount++;
	    		}
	    	}
	    	else{
	    		LaneSampleNode storedSampleNode=newAnalysisIDs.get(i).head;
	    		LaneSampleNode newSampleNode=new LaneSampleNode(flowcell, lane, sampleID,group, index, platform,"","",0,0,"","","",1,coverage,null);
	    		do
	    		{
	    			if (storedSampleNode.compareTo(newSampleNode)==true){
	    				writerErr.append("Error: "+line+" is duplicate with other records within the same analysis ID!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	    			}
	    			if (!storedSampleNode.group.equals(group)){
	    				writerErr.append("Error: "+line+" group name is different to other lane samples within one analysis ID!");
	    				writerErr.newLine();
	    				foundError=true;
	    				break;
	    			}
	    			if (!storedSampleNode.sampleID.equals(sampleID)){
	    				writerErr.append("Warning: "+line+" CGF ID is different to other lane samples within one analysis ID!");
	    				writerErr.newLine();
	    			}
	    			storedSampleNode=storedSampleNode.next;
	    		}
	    		while (storedSampleNode!=null);
	    		if (!foundError){
		    		newAnalysisIDs.get(i).add(newSampleNode);
		    		nodeCount++;
	    		}
	    	}
	    }
	    writerManifest.flush();
	    writerManifest.close();
	    System.out.println("Verification is done! Total: "+nodeCount+" lines was finished processing!");
	    reader.close();
	    writerErr.append("# Manifest verification is done! Total: "+nodeCount+" lines");
		writerErr.newLine();
		writerErr.append("");
		writerErr.newLine();
	    if (!foundError){
	    	if (!redo.equals("redo")) {
		    	System.out.println("Loading all existing merging records ...");
		    	writerErr.append("# Loading all existing merging records ...");
				writerErr.newLine();
		    	loadAllRecords(writerErr);
		    	System.out.println("Loading all existing quality trimming records ...");
		    	writerErr.append("# Loading all existing quality trimming records ...");
		    	FastqFileList fastqFileList=new FastqFileList(fileNewQualityTrimming);
		    	fastqFileList.loadFastqRecords(writerErr);
		    	
			    System.out.println("Comparing to all merging & quality trimming records ...");
			    checkAllRecords(fileListUpdate,fastqFileList);
	    	}
	    	System.out.println("Writing the output files ...");
	    	//To generate three files here: 1.generate the script for next step; 2.need to recover; 3.generate all the new or merged samples; 
	    	writeScript(fileOutput, fileListUpdate, fileNeedRestore,requestedDepth);		
	    }
	    else
	    	System.out.println("Error is found in the new Manifest!");
	    writerErr.flush();
	    writerErr.close();
	    System.out.println("Done!");

	}
	
	public boolean loadAllRecords(BufferedWriter bwErr) throws IOException{
		//AnalysisID	Flowcell	Group	Index	Platform	Lane	sample ID	Input BAM file	Modified Date	Size	Downsampled BAM file	Downsample Ratio	Log File Name
		BufferedReader readerRecords=new BufferedReader(new FileReader(fileNameCurrentAllRecords));
		String line;
		String lastAnalysisID="";
		allExistingAnalysisIDs = new ArrayList<LaneSampleList>();
		LaneSampleList analysisIDUnit=null;
		readerRecords.readLine();
		boolean foundErr=false;
		while((line=readerRecords.readLine())!=null){
			String[] fields=line.split("\t");
			//(String flowcell, String lane, String sampleID,
			//		String group, String index, String instrumentID, String filePath,String downsampledFilePath,float originalSize,
			//      float downsampledSize,String fileModifiedDate,
			//		String downsampledModifiedDate,String logFileName,float downsampleRatio,LaneSampleNode next
			
			String analysisID=fields[0];
			String[] aa=fields[1].split("_");
			String flowcell="";
			if(aa.length!=4){
				bwErr.append("Warning: the flowcell name does not contain four parts! Record line:"+line);
				bwErr.newLine();
				foundErr=true;
			    continue;
			}
			else
			  flowcell=aa[3];
			String group=fields[2];
			String index=fields[3];
			String platform=fields[4];
			String lane=fields[5];
			String sampleID=fields[6];
			String fileNameBAM=fields[7];
			String modifiedDate=fields[8];
			String fileSizeBAM=fields[9];
			String fileNameLog=fields[12];
			float downsampleRatio=Float.parseFloat(fields[11]);
			LaneSampleNode sampleNode=new LaneSampleNode(flowcell, lane, sampleID,
					group,index, platform,fileNameBAM,"",Float.parseFloat(fileSizeBAM),
					0,modifiedDate, 
					"",fileNameLog,downsampleRatio,0,null); 
			if (!analysisID.equals(lastAnalysisID)){
				if (analysisIDUnit!=null) allExistingAnalysisIDs.add(analysisIDUnit);
    		    analysisIDUnit=new LaneSampleList(analysisID,sampleNode,group);
    		    lastAnalysisID=analysisID;
			}
			else{
    		    if (analysisIDUnit.add2(sampleNode)==-1){
    		       System.err.println("Error:"+analysisID+" has duplicate lane-level BAMs!");	
    		       foundErr=true;
    		    };
    		    
			}
			
		};
		if (!foundErr){
		   if (analysisIDUnit!=null) allExistingAnalysisIDs.add(analysisIDUnit);
		}
		return foundErr;
	}
	
	
	public void checkAllRecords(String fileListUpdate, FastqFileList fastqFileList) throws IOException{
		BufferedWriter bwUpdate = new BufferedWriter(new FileWriter(
				fileListUpdate));
		ListIterator<LaneSampleList> iter=newAnalysisIDs.listIterator();
		
		while(iter.hasNext()){
			LaneSampleList currentSample=iter.next();
			boolean analysisIDFound=false;
			boolean mergeChanged=false;
			System.out.println(currentSample.analysisID);
			if (currentSample.analysisID.contains("NPC_5003_501_A"))
				System.out.println("found");
			for (int j=0;j<allExistingAnalysisIDs.size();j++){
				String currentAnalysisId=currentSample.head.group+"_"+currentSample.analysisID;
				if (currentAnalysisId.equals(allExistingAnalysisIDs.get(j).analysisID)){
					analysisIDFound=true;
			    	if (currentSample.compareTo(allExistingAnalysisIDs.get(j))){
			    		int ret=currentSample.compareToFastqFileList(fastqFileList,bwUpdate);
			    		if (ret==0){
			    			bwUpdate.append(currentSample.analysisID+" no merge needed and newly quality trimmed!");
		    			    bwUpdate.newLine();
			    			iter.remove();  // this analysisID won't do secondary analysis
			    		}
			    		else{
			    		   if (ret==1){		
			    				bwUpdate.append(currentSample.analysisID+" not all lane have been newly quality trimmed!");
			    			    bwUpdate.newLine();
			    		   }
			    			else{
			    				if (ret==2){
			    				bwUpdate.append(currentSample.analysisID+" BAM is not existed!");
			    			    bwUpdate.newLine();
			    				}
			    				else{
			    					bwUpdate.append(currentSample.analysisID+" BAM is not in new QT list!");
				    			    bwUpdate.newLine();
			    					
			    				}
			    			}
			    		}			    		
			    	}
			    	else 
			    		mergeChanged=true;
			    	break;
			    }
			}
			if (!analysisIDFound){
				bwUpdate.append(currentSample.analysisID+" is not found!");
				bwUpdate.newLine();
			}
			else
				if (mergeChanged){
				   bwUpdate.append(currentSample.analysisID+" merge changed!");
				   bwUpdate.newLine();
				}
		}
		bwUpdate.flush();
		bwUpdate.close();
	}

	public void writeScript(String fileOutput, String fileListUpdate,
			String fileNeedRestore,String requestedDepth) throws IOException, InterruptedException {
		BufferedWriter bwOutput = new BufferedWriter(new FileWriter(fileOutput));
		BufferedWriter bwUpdate = null;
		if (this.redo.equals("redo"))
		   bwUpdate = new BufferedWriter(new FileWriter(fileListUpdate));
		else
			bwUpdate= new BufferedWriter(new FileWriter(fileListUpdate,true));  //append this file
		BufferedWriter bwRestore = new BufferedWriter(new FileWriter(fileNeedRestore));
		// String suffixBAM="_HQ_paired_dedup_properly_paired_nophix";

		// if (this.dedup.equals("dedup"))
		//	suffixBAM="_dedup";
		bwOutput.append("#!/bin/sh");
		bwOutput.newLine();
		bwOutput.append("SCRIPT=$(readlink -f \"$0\")");
		bwOutput.newLine();
		bwOutput.append("SCRIPT_HOME=$(dirname \"$SCRIPT\")");
		bwOutput.newLine();
		bwOutput.append(". ${SCRIPT_HOME}/global_config_bash.rc");
		bwOutput.newLine();
		// if (this.dedup.equals("dedup")){
		 //  bwOutput.append("if [ ! -d \"$GATK_LOG_BUILD_BAM_DIR_DEDUP\" ]; then");
		//   bwOutput.newLine();
		//   bwOutput.append("mkdir $GATK_LOG_BUILD_BAM_DIR_DEDUP");
		//   bwOutput.newLine();
		//   bwOutput.append("fi");
		// }
		// else{
			bwOutput.append("if [ ! -d \"$GATK_LOG_BUILD_BAM_DIR\" ]; then");
			bwOutput.newLine();
			bwOutput.append("mkdir $GATK_LOG_BUILD_BAM_DIR");
			bwOutput.newLine();
			bwOutput.append("fi");
		// }
		bwOutput.newLine();
		
		for (int i = 0; i < newAnalysisIDs.size(); i++) {
			String downsampleRatio="1";
			String groupID=newAnalysisIDs.get(i).head.group;
			String analysisID = newAnalysisIDs.get(i).analysisID;
			LaneSampleNode sampleNodeofAnalysisID = newAnalysisIDs.get(i).head;
			
			if (newAnalysisIDs.get(i).calculateDownSampleRatios(Integer.parseInt(requestedDepth))){
				bwUpdate.append("Error: "+newAnalysisIDs.get(i).analysisID+" have two or more NextSeq!");
			    bwUpdate.newLine();
		    }
			else{
				String allBAMs = "";
				do {
					String BAM = "";
					if (sampleNodeofAnalysisID.instrumentID.toUpperCase().equals(
							"NEXTSEQ")) {
						BAM = "/DCEG/CGF/Sequencing/Illumina/NextSeq/PostRun_Analysis/Data";
						if (sampleNodeofAnalysisID.downsampleRatio>0)
						    downsampleRatio=String.valueOf(sampleNodeofAnalysisID.downsampleRatio);
						else{
							BAM="";
							bwUpdate.append("Warning: "+newAnalysisIDs.get(i).analysisID+":"+ sampleNodeofAnalysisID.instrumentID+":"+sampleNodeofAnalysisID.downsampleRatio);
							bwUpdate.newLine();
						}
					}
					else{
						if (sampleNodeofAnalysisID.instrumentID.toUpperCase().equals(
								"MISEQ")) 
				     		BAM = "/DCEG/CGF/Sequencing/Illumina/MiSeq/PostRun_Analysis/Data";
						else
							BAM = "/DCEG/CGF/Sequencing/Illumina/HiSeq/PostRun_Analysis/Data";
					}
					if (BAM.length()>0){
						BAM = BAM + "/*" + sampleNodeofAnalysisID.flowcell + "/BAM/"
								+ sampleNodeofAnalysisID.sampleID + "_"
								+ sampleNodeofAnalysisID.index + "_L00"
								+ sampleNodeofAnalysisID.lane +".bam";
						
						Process p;
					    String realName="";
					    String s="";
				        System.out.println("ls "+BAM);
				        String lscmd="ls "+BAM;
				        p = Runtime.getRuntime().exec(new String[]{"bash", "-c", lscmd});
				        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
				        while ((s = stdInput.readLine()) != null) {
				        //	System.out.println("S:"+s);
				            realName=s;
				        }
				        p.waitFor();
				        // System.out.println ("exit: " + p.exitValue());
				        p.destroy();
				        // System.out.println(realName);
				        File bamFile=new File(realName.trim());
				        Path realPath = bamFile.toPath().toRealPath();
				        if (realPath.toString().toUpperCase().contains("DOWNSAMPLE")) {
				        	bwUpdate.append("Warning: "+realPath.toString()+" has been downsampled!");
							bwUpdate.newLine();
				        	
				        }
				        	
				        		
						if (realName.trim().length()==0)
						{
							// System.out.println ("Not found!");
							bwRestore.append(BAM);
							bwRestore.newLine();	
							allBAMs = allBAMs + BAM + " ";
						}
						else{
							// System.out.println ("Found!");
						    allBAMs = allBAMs + realName + " ";
						}
					}
					sampleNodeofAnalysisID = sampleNodeofAnalysisID.next;
				} 
				while (sampleNodeofAnalysisID != null);
				
				bwUpdate.append(groupID+"_"+analysisID + "\t" + allBAMs);
				bwUpdate.newLine();
				if (downsampleRatio.equals("1.0"))  downsampleRatio="1";
				// if (dedup.equals("dedup"))
				//	bwOutput.append("CMD=\"qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR_DEDUP}/_build_bam_"+groupID+"_"+analysisID+".stdout -e ${GATK_LOG_BUILD_BAM_DIR_DEDUP}/_build_bam_"+groupID+"_"+analysisID+".stderr "
		        //            +"-N BAMCOMB."+groupID+"_"+analysisID+" -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v3.sh "+groupID+"_"+analysisID +" "+variantType+" "+downsampleRatio+" " + allBAMs+"\"");
				// else
				System.out.println("v4"+this.fileName);
				bwOutput.append("CMD=\"qsub -q $QUEUE -o ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_"+groupID+"_"+analysisID+".stdout -e ${GATK_LOG_BUILD_BAM_DIR}/_build_bam_"+groupID+"_"+analysisID+".stderr "
                    +"-N BAMCOMB."+groupID+"_"+analysisID+" -v SCRIPT_HOME=$SCRIPT_HOME -S /bin/sh ${SCRIPT_HOME}/gatk_build_bam_for_single_name_v4.sh "+groupID+"_"+analysisID +" "+downsampleRatio+" " + this.fileName+" "+ allBAMs+"\"");
				bwOutput.newLine();
				bwOutput.append("echo $CMD");
				bwOutput.newLine();
				bwOutput.append("eval $CMD");
				bwOutput.newLine();
			}

		}
		bwUpdate.flush();
		bwUpdate.close();
		bwOutput.flush();
		bwOutput.close();
		bwRestore.flush();
		bwRestore.close();
		System.out.println("Done!");
	}
	
}


