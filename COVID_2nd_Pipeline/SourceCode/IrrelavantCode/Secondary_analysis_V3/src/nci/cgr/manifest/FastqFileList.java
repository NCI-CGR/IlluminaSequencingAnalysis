package nci.cgr.manifest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class FastqFileList {

	public FastqFileList(String filename) {
		super();
		this.filename = filename;
		List<FastqFile> fastqFileList=new ArrayList<FastqFile>();
		this.fastqFiles=fastqFileList;
	}
	
	
	public int loadFastqRecords(BufferedWriter bwErr) throws IOException, ParseException{
		@SuppressWarnings("resource")
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		String line="";
		int lineNum=1;
		int ret=0;
	    while ((line=reader.readLine())!=null){
	    	String aa[]=line.split("\\t");
	    	if (aa[3].trim().length()==0 ||  !aa[1].contains("_") 
	    			|| !aa[5].contains("_") || aa[2].trim().length()==0 || aa[4].trim().length()==0  
	    			|| aa[0].trim().length()==0){
	    		bwErr.append("Error: FASTQ File list format is error on line: " + lineNum+" ("+line+")");
	    		bwErr.newLine();
	    		bwErr.flush();
	    		ret=-1;
	    	}
	    	else{
	    		if (aa[6].startsWith("NEW")) {
		    	//	System.out.println(aa[0]+"\t"+aa[1]+"\t"+aa[2]+"\t"+aa[3]+"\t"+aa[4]+"\t"+aa[5]);
			    	String lane=aa[3].substring(1, 2);
			    
			    	String fields[]=aa[1].split("_");
			    	
			    	String sampleFields[]=aa[5].split("_");
			    	
			    	if (fields.length!=4 || sampleFields.length<2) {
			    		bwErr.append("Error: FASTQ File list format is error on line: " + lineNum+" ("+line+")");
		    	    	bwErr.newLine();
		    	    	bwErr.flush();
		    	    	ret=-1;
			    	}
			    	else{
			    		String expectedPattern="yyyy-MM-dd";
			    		SimpleDateFormat formatter=new SimpleDateFormat(expectedPattern);
			    		Date date=formatter.parse(aa[8]);
			    		String sampleName=sampleFields[1];
			    		if (sampleFields.length>2){
			    		   for (int i=2;i<sampleFields.length;i++)
			    		      sampleName=sampleName+"_"+sampleFields[i];
			    		 }
			    
				    	FastqFile fastqFile=new FastqFile(aa[0],fields[3],lane,aa[2],aa[4],sampleName,date);
				    	fastqFiles.add(fastqFile);
			    	}
	    		}
	    	}
			lineNum++;
		}
	    return ret;
	}
	
	String filename;
	List<FastqFile> fastqFiles;
	
}