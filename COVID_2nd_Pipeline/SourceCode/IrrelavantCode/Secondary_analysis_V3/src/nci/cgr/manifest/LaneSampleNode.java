package nci.cgr.manifest;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Date;

public class LaneSampleNode {
	public LaneSampleNode(String flowcell, String lane, String sampleID,
			String group, String index, String instrumentID, String filePath,String downsampledFilePath,float originalSize,float downsampledSize,String fileModifiedDate,
			String downsampledModifiedDate,String logFileName,float downsampleRatio,float coverage,LaneSampleNode next) {
		super();
		this.flowcell = flowcell;
		this.lane = lane;
		this.sampleID = sampleID;
		this.group = group;
		this.index = index;
		this.instrumentID=instrumentID;
		this.filePath=filePath;
		this.downsampledFilePath=downsampledFilePath;
		this.fileSize=originalSize;
		this.downsampledFileSize=downsampledSize;
		this.downsampledModifiedDate=downsampledModifiedDate;
		this.fileModifiedDate=fileModifiedDate;
		this.logFileName=logFileName;
		this.downsampleRatio=downsampleRatio;
		this.coverage=coverage;
		this.next = null;
	}
	
	
	public LaneSampleNode(String flowcell, String lane, String sampleID,
			String group, String index, String instrumentID, String filePath,String downsampledFilePath,float originalSize,float downsampledSize,String fileModifiedDate,
			String downsampledModifiedDate,String logFileName,float downsampleRatio,float coverage,LaneSampleNode next,Date seqDate) {
		super();
		this.flowcell = flowcell;
		this.lane = lane;
		this.sampleID = sampleID;
		this.group = group;
		this.index = index;
		this.instrumentID=instrumentID;
		this.filePath=filePath;
		this.downsampledFilePath=downsampledFilePath;
		this.fileSize=originalSize;
		this.downsampledFileSize=downsampledSize;
		this.downsampledModifiedDate=downsampledModifiedDate;
		this.fileModifiedDate=fileModifiedDate;
		this.logFileName=logFileName;
		this.downsampleRatio=downsampleRatio;
		this.coverage=coverage;
		this.seqDate=seqDate;
		this.next = null;
	}
   
	String flowcell;
	String lane;
	String sampleID;
	String group;
	String index;
	String instrumentID;
	String filePath;
	String downsampledFilePath;
	float fileSize;
	float downsampledFileSize;
	String fileModifiedDate;
	String downsampledModifiedDate;
	String logFileName;
	float downsampleRatio;
	float coverage;
	Date seqDate;
	LaneSampleNode next;
	
	
	
	
	
	public String getFlowcell() {
		return flowcell;
	}
	public void setFlowcell(String flowcell) {
		this.flowcell = flowcell;
	}
	public String getLane() {
		return lane;
	}
	public void setLane(String lane) {
		this.lane = lane;
	}
	public String getSampleID() {
		return sampleID;
	}
	public void setSampleID(String sampleID) {
		this.sampleID = sampleID;
	}
	public String getGroup() {
		return group;
	}
	public void setGroup(String group) {
		this.group = group;
	}
	public String getIndex() {
		return index;
	}
	public void setIndex(String index) {
		this.index = index;
	}
	public String getInstrumentID() {
		return instrumentID;
	}
	public void setInstrumentID(String instrumentID) {
		this.instrumentID = instrumentID;
	}
	
	public boolean compareTo(LaneSampleNode other){
		if ( this.index.equals(other.index) && this.lane.equals(other.lane)
				&& this.flowcell.equals(other.flowcell)&& this.sampleID.equals(other.sampleID))
			return true;
		else
			return false;
	}
	
	public boolean compareTo3(LaneSampleNode other){
		if ( this.index.equals(other.index) && this.lane.equals(other.lane)
				&& ((other.flowcell.length()==this.flowcell.length()+1) && (other.flowcell.contains(this.flowcell))) ||  other.flowcell.equals(this.flowcell))
//			if ( this.index.equals(other.index) && this.lane.equals(other.lane))
//					&& ((other.flowcell.length()==this.flowcell.length()+1) && (other.flowcell.contains(this.flowcell))) ||  other.flowcell.equals(this.flowcell))	
		return true;
		else
			return false;
	}
	
	public boolean compareTo2(LaneSampleNode other){
		if (  this.lane.equals(other.lane)){
			System.out.println(this.lane+" "+other.lane);
			if ((other.flowcell.length()==this.flowcell.length()+1) && (other.flowcell.contains(this.flowcell))){
			 // System.out.println(this.flowcell+" "+other.flowcell);  
			  return true;
			}
			else
				return false;
		}
		else{
			if(other.instrumentID.equals("NextSeq")){
				if ((other.flowcell.length()==this.flowcell.length()+1) && (other.flowcell.contains(this.flowcell)) && 
						(this.lane.equals("1") || this.lane.equals("2") || this.lane.equals("3") || this.lane.equals("4")))
					return true;
				else
					return false;
			}
			else
			  return false;
		}
	}
	
	public boolean compareToFastq(FastqFile fastqFile,BufferedWriter bw) throws IOException {
		// bw.append("FASTQ:"+fastqFile.index+" "+fastqFile.lane+" "+fastqFile.flowcell+" "+fastqFile.sample+" "+fastqFile.platform);
		// bw.append("BAM:"+this.index+" "+this.lane+" "+this.flowcell+" "+this.sampleID+" "+this.instrumentID);
		// bw.newLine();
		// bw.flush();
		// if (this.sampleID.equals(fastqFile.sample))
		//	System.out.println("found!");
		if (this.index.equals(fastqFile.index) && this.lane.equals(fastqFile.lane) && this.flowcell.equals(fastqFile.flowcell) && this.sampleID.equals(fastqFile.sample) && this.instrumentID.equals(fastqFile.platform))
			  return true;
		else
		      return false;
		
	}
	
	public String getFilePath() {
		return filePath;
	}
	public void setFilePath(String filePath) {
		this.filePath = filePath;
	}
	public String getDownsampledFilePath() {
		return downsampledFilePath;
	}
	public void setDownsampledFilePath(String downsampledFilePath) {
		this.downsampledFilePath = downsampledFilePath;
	}
	public float getFileSize() {
		return fileSize;
	}
	public void setFileSize(float fileSize) {
		this.fileSize = fileSize;
	}
	public float getDownsampledFileSize() {
		return downsampledFileSize;
	}
	public void setDownsampledFileSize(float downsampledFileSize) {
		this.downsampledFileSize = downsampledFileSize;
	}
	public String getFileModifiedDate() {
		return fileModifiedDate;
	}
	public void setFileModifiedDate(String fileModifiedDate) {
		this.fileModifiedDate = fileModifiedDate;
	}
	public String getDownsampledModifiedDate() {
		return downsampledModifiedDate;
	}
	public void setDownsampledModifiedDate(String downsampledModifiedDate) {
		this.downsampledModifiedDate = downsampledModifiedDate;
	}
	public String getLogFileName() {
		return logFileName;
	}
	public void setLogFileName(String logFileName) {
		this.logFileName = logFileName;
	}
	public float getDownsampleRatio() {
		return downsampleRatio;
	}
	public void setDownsampleRatio(float downsampleRatio) {
		this.downsampleRatio = downsampleRatio;
	}


	public Date getSeqDate() {
		return seqDate;
	}


	public void setSeqDate(Date seqDate) {
		this.seqDate = seqDate;
	}
	
	
}
