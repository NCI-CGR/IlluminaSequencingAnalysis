package nci.cgr.manifest;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.*;
import static java.nio.file.FileVisitResult.*;
import static java.nio.file.FileVisitOption.*;

public class LaneSampleList {
	public LaneSampleList(String analysisID, LaneSampleNode laneSampleNode, String group) {
		super();
		this.analysisID = analysisID;
		this.group=group;
		this.head = laneSampleNode;
		String bamDir="/DCEG/Projects/Exome/SequencingData/BAM_reformatted/BAM_original";
		File file=new File(bamDir+"/"+group+"/"+group+"_"+analysisID+".bam");
		if (file.exists())
		   this.bamUpdateDate=new Date(file.lastModified());
		else
			this.bamUpdateDate=null;
	}
	String analysisID;
	String group;
	LaneSampleNode head;
	Date bamUpdateDate;
	
	public void add(LaneSampleNode newNode){
		LaneSampleNode currentNode=head;
		if (currentNode==null)
			this.head=newNode;
		else{
			while(currentNode.next!=null){
				currentNode=currentNode.next;
			}
			currentNode.next=newNode;
		}
	}
	
	public int add2(LaneSampleNode newNode){
		LaneSampleNode currentNode=head;
		if (currentNode==null){
			this.head=newNode;
		    return 1;	
		}
		else{
			while(currentNode.next!=null){
				if (newNode.compareTo(currentNode)==true){
					return -1;
				}
				else
				  currentNode=currentNode.next;
			}
			currentNode.next=newNode;
			return 1;
		}
	}
	
	public int length(){
		LaneSampleNode node=head;
		int length=0;
		while(node!=null){
			length++;
			node=node.next;
		}
		return length;
	}
	
	public boolean calculateDownSampleRatios(float requiredDepth) {
		// return true if two or more NextSeq is found else return false
		LaneSampleNode node = head;
		float totalCoverageNoNextSeq = 0;
		boolean foundNextSeq = false;
		
		int nextSeqCount=0;
		while (node != null) {
			if (node.instrumentID.toUpperCase().equals("NEXTSEQ")){
				foundNextSeq = true;
				nextSeqCount++;
			}
			else
				totalCoverageNoNextSeq += node.coverage;
			node = node.next;
		}
		
		if (foundNextSeq) {
			if (nextSeqCount>1)
				return true;
			else{
				if (totalCoverageNoNextSeq > requiredDepth) {
					node = head;
					while (node != null) {
						if (node.instrumentID.toUpperCase().equals("NEXTSEQ"))
							node.downsampleRatio = 0;
						node = node.next;
					}
				} else {
					if (totalCoverageNoNextSeq > requiredDepth / 2) {
						node = head;
						while (node != null) {							
							if (node.instrumentID.toUpperCase().equals("NEXTSEQ")) {
								if (node.coverage > requiredDepth
										- totalCoverageNoNextSeq)
									node.downsampleRatio = (requiredDepth - totalCoverageNoNextSeq)	/ node.coverage;
								else
									node.downsampleRatio = 1;
							}
							node = node.next;
						}
					} else {
						node = head;
						while (node != null) {
							if (node.instrumentID.toUpperCase().equals("NEXTSEQ"))
								if (totalCoverageNoNextSeq>0)
									if (totalCoverageNoNextSeq<node.coverage)
								       node.downsampleRatio = totalCoverageNoNextSeq / node.coverage;
									else
										node.downsampleRatio=1;
								else
								   node.downsampleRatio = 1;
							node = node.next;							
						}	
					}
				}
				return false;
			}
		}
		else
			return false;
	}
	
	
	public boolean compareTo(LaneSampleList anotherOne) {
		LaneSampleNode node1 = this.head;
		LaneSampleNode node2 = null;
		int length1 = this.length();
		int length2 = anotherOne.length();
		if (length1 != length2)
			return false;
		else {
			int cmpLength = 0;
			do {
				node2 = anotherOne.head;
				Boolean found=false;
				do {
					if (node1.compareTo(node2) == true) {
						found=true;
						break;
					}
					node2=node2.next;
				} while (node2 != null);
				if (found) cmpLength++;
				else
					return false;
				node1 = node1.next;
			} while (node1 != null);
			if (cmpLength == length1)
				return true;
			else
				return false;

		}
	}
	
	public boolean compareTo3(LaneSampleList anotherOne) {
		LaneSampleNode node1 = this.head;
		LaneSampleNode node2 = null;
		int length1 = this.length();
		int length2 = anotherOne.length();
		if (length1 != length2)
			return false;
		else {
			int cmpLength = 0;
			do {
				node2 = anotherOne.head;
				Boolean found=false;
				do {
					if (node1.compareTo3(node2) == true) {
						found=true;
						break;
					}
					node2=node2.next;
				} while (node2 != null);
				if (found) cmpLength++;
				else
					return false;
				node1 = node1.next;
			} while (node1 != null);
			if (cmpLength == length1)
				return true;
			else
				return false;

		}
	}
	
	public int compareToFastqFileList(FastqFileList fastqFileList,BufferedWriter bw) throws IOException{
		// Return 1 if modified date of any FASTQ is after BAM update date or BAM is not existed
		// Return 2 if no BAM existed;
		// Return 0 if BAM is newer than all underlying FASTQs;
		// Return 3 if BAM has some FASTQ is not in new trimming list;
		
		Date bamUpdateDate=this.bamUpdateDate;
		int ret=0;
		if (bamUpdateDate!=null){
			LaneSampleNode node=this.head;
			while(node!=null){
				boolean found=false;
				for (int i=0;i<fastqFileList.fastqFiles.size();i++){
				  if(node.compareToFastq(fastqFileList.fastqFiles.get(i),bw)){
					  found=true;
					  if (fastqFileList.fastqFiles.get(i).updateDate.compareTo(bamUpdateDate)>0){ //FASTQ newer than BAM
						  ret=1;
						  DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
			               
			                //to convert Date to String, use format method of SimpleDateFormat class.
			              String bamDate = dateFormat.format(this.bamUpdateDate);
			              String fastqDate=dateFormat.format(fastqFileList.fastqFiles.get(i).updateDate);
						  bw.append("Found:"+this.group+"/"+this.analysisID+" date:"+bamDate+	" "+fastqFileList.fastqFiles.get(i).flowcell+"/"+
			            		  fastqFileList.fastqFiles.get(i).index+"/"+fastqFileList.fastqFiles.get(i).lane+"/"+fastqFileList.fastqFiles.get(i).sample
			            		  +" date:"+fastqDate);
			              bw.newLine();
			              bw.flush();
					   
					  }
					  break;
				  }
				}
				if (!found){
					  ret=3;
					  bw.append(node.flowcell+"/"+node.index+"/"+node.sampleID+"/"+node.lane+" is not in new QT list!");
					  bw.newLine();
			    }
				if ((ret==1) || (ret==3)) break;
				node=node.next;
			}
		}
		else{
			ret=2;
			bw.append(this.group+"/"+this.analysisID+" (BAM) is not existed!");
			bw.newLine();
			bw.flush();
		}
		return ret;
	}
	
	
	public boolean compareToFastqFileList(FastqFileList fastqFileList, BufferedWriter bwOutput,BufferedWriter bwRedo) throws IOException, ParseException, InterruptedException{
    // hard coded two specific cases for printout the fastq list		
    // This function is for comparing the manifest and new QT list to determine which FASTQ needed restoration 
		LaneSampleNode node1 = this.head;
		String seqDir="/DCEG/CGF/Sequencing/Illumina/";
		
		boolean allLanesNewTrimmed=true;
		do{
			boolean found=false;
			for (int i=0;i<fastqFileList.fastqFiles.size();i++){
				if (node1.compareToFastq(fastqFileList.fastqFiles.get(i),bwOutput)){
					found=true;
				    break;
				}
			}
			if (!found){
				allLanesNewTrimmed=false;
		    	String expectedPattern = "MM/dd/yyyy";
		        SimpleDateFormat formatter = new SimpleDateFormat(expectedPattern);
		        Date thresholdDate=formatter.parse("04/30/2014");  	    	
				String fileNames="";
				boolean notExists=false;
		    	if (node1.getSeqDate().before(thresholdDate)){
		    		if (node1.getFlowcell().equals("BC3CLRACXX"))
		    			fileNames=seqDir+node1.getInstrumentID().trim()+"/PostRun_Analysis/Data/*"+node1.getFlowcell().trim()+"/CASAVA_1mismatch/Project*/Sample_"+node1.getSampleID().trim()+"/"+node1.getSampleID().trim()+"_"+node1.getIndex().trim()+"_L00"+node1.getLane().trim()+"_R*_001.fastq.gz";
		    		else
		    			
		    			fileNames=seqDir+node1.getInstrumentID().trim()+"/PostRun_Analysis/Data/*"+node1.getFlowcell().trim()+"/CASAVA/Project*/Sample_"+node1.getSampleID().trim()+"/"+node1.getSampleID().trim()+"_"+node1.getIndex().trim()+"_L00"+node1.getLane().trim()+"_R*_001.fastq.gz";	
		    	}
		    	else{
		    		if (node1.getFlowcell().equals("AC6PVCANXX"))
		    			fileNames=seqDir+node1.getInstrumentID().trim()+"/PostRun_Analysis/Data/*"+node1.getFlowcell().trim()+"/CASAVA/L"+node1.getLane().trim()+"/Project*/Sample_"+node1.getSampleID().trim()+"/"+node1.getSampleID().trim()+"__L00"+node1.getLane().trim()+"_R*_001.fastq.gz";
		    		else
		    		    fileNames=seqDir+node1.getInstrumentID().trim()+"/PostRun_Analysis/Data/*"+node1.getFlowcell().trim()+"/CASAVA/L"+node1.getLane().trim()+"/Project*/Sample_"+node1.getSampleID().trim()+"/"+node1.getSampleID().trim()+"_"+node1.getIndex().trim()+"_L00"+node1.getLane().trim()+"_R*_001.fastq.gz";
		    	}
		    	Process p;
		    	
    		    // p = Runtime.getRuntime().exec("ls "+fileNames.trim());
    		    p = Runtime.getRuntime().exec(new String[] { "/bin/sh", "-c", "ls "+fileNames.trim()});
    		    p.waitFor();
    		    BufferedReader buf = new BufferedReader(new InputStreamReader(p.getErrorStream()));
    		    
    		    String line = "";
    		    while ((line = buf.readLine()) != null) {
    		      System.out.println("return Error:"+line);
    		      if (line.contains("cannot")){
    		    	  notExists=true;		    		    	  
    		      }
    		    }
    		    
    		    BufferedReader bufInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
    		    line="";
    		    while ((line = bufInput.readLine()) != null) {
        		      System.out.println("return Input:"+line);
        		}
    		    
    		    System.out.println("ls "+fileNames.trim());
    		    
    		    
    		    
		    	if (node1.getSeqDate().before(thresholdDate)  || node1.getSeqDate().equals(thresholdDate)){
		    	    if (node1.getFlowcell().equals("BC3CLRACXX")){
		    		    if (notExists){
		    		    	bwOutput.append(fileNames);
				     		bwOutput.newLine();
			    	    	bwOutput.flush();
		    		    }
			        	bwRedo.append(seqDir+node1.getInstrumentID()+"/PostRun_Analysis/Data/*"+node1.getFlowcell()+"/CASAVA_1mismatch/Project*/Sample_"+node1.getSampleID());
			    	    bwRedo.newLine();
			    	    bwRedo.flush(); 
		    		    
		    		}
		    		else{
		    		    if (notExists){
		    		    	bwOutput.append(fileNames);
				     		bwOutput.newLine();
			    	    	bwOutput.flush();
		    		    }
			    	    bwRedo.append(seqDir+node1.getInstrumentID()+"/PostRun_Analysis/Data/*"+node1.getFlowcell()+"/CASAVA/Project*/Sample_"+node1.getSampleID());
			    	    bwRedo.newLine();
			    	    bwRedo.flush(); 
		    		    				    	
		    		}
		    	}
			    else{
			        if (notExists){
				      bwOutput.append(fileNames);
				      bwOutput.newLine();
	    	     	  bwOutput.flush();
	    		    }
	    	     	bwRedo.append(seqDir+node1.getInstrumentID()+"/PostRun_Analysis/Data/*"+node1.getFlowcell()+"/CASAVA/L"+node1.getLane()+"/Project*/Sample_"+node1.getSampleID());
	    	    	bwRedo.newLine();
	    	    	bwRedo.flush();
	    		    
			    }
			}
	    	node1=node1.next;
		}while (node1 !=null);
		return allLanesNewTrimmed;
	}
	
	public boolean compareTo2(LaneSampleList anotherOne) {
		LaneSampleNode node1 = this.head;
		LaneSampleNode node2 = null;
		int length1 = this.length();
		int length2 = anotherOne.length();
		if (length1 != length2){
			
			LaneSampleNode node3 = anotherOne.head;
			Boolean found=false;
			do {
				if (node3.instrumentID.equals("NextSeq")) {
					found=true;
					break;
				}
				node3=node3.next;
			} while (node3 != null);
			if (found){
				int cmpLength = 0;
				do {
					node2 = anotherOne.head;
					found=false;
					do {
						if (node1.compareTo2(node2) == true) {
							found=true;
						//	System.out.println("found!");
							break;
						}
						node2=node2.next;
					} while (node2 != null);
					if (found) cmpLength++;
					else
						return false;
					node1 = node1.next;
				} while (node1 != null);
				if (cmpLength == length1){
					// System.out.println("found Match!");
					return true;
				}
				else{
					// System.out.println("Length different");
					return false;
				}
			}
				
			else
				return false;
		}
		else {
			int cmpLength = 0;
			do {
				node2 = anotherOne.head;
				Boolean found=false;
				do {
					if (node1.compareTo2(node2) == true) {
						found=true;
						// System.out.println("found!");
						break;
					}
					node2=node2.next;
				} while (node2 != null);
				if (found) cmpLength++;
				else
					return false;
				node1 = node1.next;
			} while (node1 != null);
			if (cmpLength == length1){
				// System.out.println("found Match!");
				return true;
			}
			else{
				// System.out.println("Length different");
				return false;
			}

		}
	}
}