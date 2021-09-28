package nci.cgr.manifest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

public class BamInfo {
	public BamInfo() {
		allExistingBamAnalysisIDs=null;
	}

	List<LaneSampleList> allExistingBamAnalysisIDs;
	public void parseBamInfo(String fileRecords,String fileBamOutput,String fileReport,String fileErr) throws IOException{
		Manifest manifest=new Manifest(fileRecords,"","");
		//BufferedWriter bwErr=new BufferedWriter(new FileWriter(fileErr));
		manifest.verifyManifest(fileErr);
		//bwErr.flush();
		//bwErr.close();
		
		// Bam parsing results
		loadAllBamRecords(fileBamOutput);
		
		
		// Compare the manifest is consistency with the mappings in the BAMs
		compareTwoRecordLists(manifest.newAnalysisIDs,fileReport);
	}
	
	public void compareTwoRecordLists(List<LaneSampleList> analysisIDs,String fileReport) throws IOException{
        ListIterator<LaneSampleList> iter=allExistingBamAnalysisIDs.listIterator(); //BAM parse
		BufferedWriter bwReport=new BufferedWriter(new FileWriter(fileReport));
		System.out.println("Comparing two lists ...");
		if (allExistingBamAnalysisIDs.size()!=analysisIDs.size()) {
			boolean found=false;
			for (int j=0;j<analysisIDs.size();j++) {
				found=false;
				// System.out.println(j);
				for (int i=0;i<allExistingBamAnalysisIDs.size();i++) {
				//	System.out.println("i"+i);
					String tmpID=analysisIDs.get(j).group+"_"+analysisIDs.get(j).analysisID;
					if (tmpID.equals(allExistingBamAnalysisIDs.get(i).analysisID)){
						found=true;
						break;
					}
				}
				if (!found) {
					bwReport.append("Error:"+analysisIDs.get(j).analysisID+" is not found in the BAM Parsing result!");
					bwReport.newLine();
				}
			}
			
		}
		while(iter.hasNext()){
			LaneSampleList currentSample=iter.next();
			boolean analysisIDFound=false;
			boolean mergeChanged=false;
			for (int j=0;j<analysisIDs.size();j++){
				String currentAnalysisId=currentSample.analysisID;
				System.out.println(currentAnalysisId);
				System.out.println(analysisIDs.get(j).analysisID);
				System.out.println(analysisIDs.get(j).group);
//				if (currentAnalysisId.equals("CHR_CHR_4124_1001_A")){						
//						System.out.println("Found !");
//			    }
				if (currentAnalysisId.equals(analysisIDs.get(j).group+"_"+analysisIDs.get(j).analysisID)){
					analysisIDFound=true;					
					
					
			    	if (currentSample.compareTo3(analysisIDs.get(j))){
			    		bwReport.append(currentSample.analysisID+" is matched!");
						bwReport.newLine();
			    		iter.remove();
			    	}
			    	else 
			    		mergeChanged=true;
			    	break;
			    }
			}
			if (!analysisIDFound){
				bwReport.append("Error:"+currentSample.analysisID+" is not found in the merge log!");
						
				bwReport.newLine();
			}
			else
				if (mergeChanged){
				   bwReport.append("Error:"+currentSample.analysisID+" merge changed!");
				   bwReport.newLine();
				}
		}
		bwReport.flush();
		bwReport.close();
		System.out.println("Comparison is done!");
	}
	
	public void loadAllBamRecords(String fileBamOutput) throws IOException{
		BufferedReader brBamOutput=new BufferedReader(new FileReader(fileBamOutput));
		String line="";
		String lastAnalysisID="";
		LaneSampleList analysisIDUnit=null;
		allExistingBamAnalysisIDs=new ArrayList<LaneSampleList>();;
	    while((line=brBamOutput.readLine())!=null){
			String[] fields=line.split("\t");	
			String analysisID=fields[0];
			String flowcell=fields[2];
			String machine=fields[1];
			String lane=fields[3];
			String index=fields[4];
			// String[] elements=pu.split("_");
			// String index=elements[elements.length-4];
			String platform="HiSeq";
			if (machine.equals("K00278"))
				platform="Hiseq";
			LaneSampleNode sampleNode=new LaneSampleNode(flowcell, lane, "",
					"",index, platform,"","",0,0,"", "","",0,0,null); 
			if (!analysisID.equals(lastAnalysisID)){
				if (analysisIDUnit!=null) allExistingBamAnalysisIDs.add(analysisIDUnit);
    		    analysisIDUnit=new LaneSampleList(analysisID,sampleNode,"NA");
    		    lastAnalysisID=analysisID;
			}
			else
    		    analysisIDUnit.add(sampleNode);
		};
		if (analysisIDUnit!=null) allExistingBamAnalysisIDs.add(analysisIDUnit);
		brBamOutput.close();
	}
}
