package nci.cgr.manifest;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;

public class LogMerging {
	public LogMerging(String logDir) {
	//	super();
		this.logDir = logDir;
	}

	
	String logDir = "";

	public String getLogDir() {
		return logDir;
	}

	public void setLogDir(String logDir) {
		this.logDir = logDir;
	}

	
	@SuppressWarnings("unchecked")
	public void scanLogs(String fileOut, String fileErr) throws IOException, InterruptedException {
		File f = new File(logDir);

		File[] files = f.listFiles();
		FileWriter fileErrHandler = new FileWriter(fileErr,false);
		BufferedWriter writerErr = new BufferedWriter(fileErrHandler);

		FileWriter fileOutHandler = new FileWriter(fileOut,false);
		BufferedWriter writerOut = new BufferedWriter(fileOutHandler);

		// Arrays.sort( files, new Comparator()
		// {
		// public int compare(Object o1, Object o2) {
		//
		// if (((File)o1).lastModified() > ((File)o2).lastModified()) {
		// return -1;
		// } else if (((File)o1).lastModified() < ((File)o2).lastModified()) {
		// return +1;
		// } else {
		// return 0;
		// }
		// }
		//
		// });
		List<LaneSampleList> analysisIDs = new ArrayList<LaneSampleList>();
		List<File> holder = new ArrayList<File>();

		for (File ff : files) {

			System.out.println(ff.getAbsolutePath());
			if (ff.isDirectory()) {
				File[] logFiles = ff.listFiles();
				for (int j = 0; j < logFiles.length; j++) {
					System.out.println(ff.getAbsolutePath() + " "
							+ logFiles[j].getAbsolutePath());
					if (logFiles[j].getAbsolutePath().endsWith("stdout")) {
						System.out.println("adding "
								+ logFiles[j].getAbsolutePath());
						holder.add(logFiles[j]);
					}
				}
			}
		}
		System.out.println("Sorting total number:" + holder.size() + " ...");

		Collections.sort(holder, new Comparator<File>() {
			public int compare(File o1, File o2) {
				return Long.compare(o1.lastModified(), o2.lastModified());
			}
		});

		for (int i = holder.size() - 1; i >= 0; i--) {
			Date d = new Date(holder.get(i).lastModified());
			System.out.println("Processing  "+i+":"+ holder.get(i).getAbsolutePath()
					+ "\t" + holder.get(i).lastModified() + "\t" + d);

			writerErr.write("Processing " + holder.get(i).getAbsolutePath()
					+ "\t" + d);
			writerErr.newLine();
			FileReader logFile = new FileReader(holder.get(i).getAbsolutePath());
			@SuppressWarnings("resource")
			BufferedReader readerLog2 = new BufferedReader(logFile);
			String line = "";
			boolean errorFound=false;
			while ((line = readerLog2.readLine()) != null) {
				if (line.contains("Error")) {
					writerErr.write("Error: "
							+ holder.get(i).getAbsolutePath()
							+ " Merging log contains some errors!");
					writerErr.newLine();
					errorFound=true;
					break;
				}				
			}
			readerLog2.close();
			logFile.close();
			System.out.println(errorFound);
			if (errorFound) continue;
			System.out.println("Starting checking"+holder.get(i).getAbsolutePath());
			BufferedReader readerLog = new BufferedReader(new FileReader(holder.get(i).getAbsolutePath()));
		
			int countCleaning = 0;
			float downsampleRatio=1;
			while ((line = readerLog.readLine()) != null) {
				if (line.contains("PROBABILITY=")){
					String[] aa=line.split("PROBABILITY=");
					downsampleRatio=Float.parseFloat(aa[1].trim());
				}
				if (line.startsWith("Building")) {
					String[] aa = line.split(" ");
					String analysisID = aa[4];

					String[] parts = analysisID.split("_");
					String group = parts[0];
					boolean found = false;
					for (int j = 0; j < analysisIDs.size(); j++) {
						if (analysisID.equals(analysisIDs.get(j).analysisID)) {
							found = true;
							break;
						}
					}
					if (found)	continue;
					
					LaneSampleList ll = new LaneSampleList(analysisID, null,group);
					boolean parsingErr = false;

					line = readerLog.readLine();
					System.out.println(line);
					if (line.startsWith("Input:")) {
						String lineDownsampled = readerLog.readLine();
						System.out.println(lineDownsampled);
						if (!lineDownsampled.startsWith("Downsampled") && !lineDownsampled.startsWith("Output")) {
								System.out.println("Error");
						    	writerErr.write("Error: "
									+ holder.get(i).getAbsolutePath()
									+ " Input is not followed by Downsampled or Output!");
						    	writerErr.newLine();
						    	parsingErr = true;
						    	// System.out.println("going break!");
						    	break;							
						} else {
							String[] parts1=line.split("Input:");
							String[] bamFiles =parts1[1].split(" ");
							for (int kk = 1; kk < bamFiles.length; kk++) {
								Path originalDownsampledPath = null;
								float downsampledSize=0;
								String downsampledFileDate="";
								if (lineDownsampled.startsWith("Downsampled")) {
									String[] parts2 = lineDownsampled.split("Input:");
									String[] downSampledFiles = parts2[1].split(" ");
									if (downSampledFiles.length != bamFiles.length) {
										writerErr.write("Error: "+ holder.get(i).getAbsolutePath()+ " the numbers of Input and downsampled input are not matched!");
										writerErr.newLine();
										parsingErr = true;
										break;
									}
									if (!downSampledFiles[kk].equals(bamFiles[kk])) {  //Now is only for NextSeq
										Path downsampledFilekk = FileSystems.getDefault().getPath(downSampledFiles[kk]);
										if (Files.exists(downsampledFilekk)) {
											if (Files.isSymbolicLink(downsampledFilekk)) 
												originalDownsampledPath = Files.readSymbolicLink(downsampledFilekk);
											else
												originalDownsampledPath = downsampledFilekk;
											if (Files.exists(originalDownsampledPath)) {
												downsampledSize = Files.size(originalDownsampledPath);
												downsampledFileDate = Files.getLastModifiedTime(originalDownsampledPath).toString();
											} else {
												String s = "";
												String realName = "";
												Process p;
												p = Runtime.getRuntime().exec("readlink -f "+ downSampledFiles[kk]);
												BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
												while ((s = br.readLine()) != null) 
													realName = s;
												
												p.waitFor();
												// System.out.println("exit: "	+ p.exitValue());
												p.destroy();

												if (Files.exists(FileSystems.getDefault().getPath(realName))) {
													originalDownsampledPath = FileSystems.getDefault().getPath(realName);
													downsampledSize = Files.size(FileSystems.getDefault().getPath(realName));
													downsampledFileDate = Files.getLastModifiedTime(FileSystems.getDefault().getPath(realName)).toString();
												} else {
													downsampledSize = 0;
													downsampledFileDate = "Original file is not exisited";
												}
											}
										} else {
											downsampledSize = 0;
											downsampledFileDate = "Not exisited";
											originalDownsampledPath = downsampledFilekk;
										}
									}
								}
								String originalFileDate="";
								Path filekk = FileSystems.getDefault().getPath(bamFiles[kk]);
								Path originalPath = null;
								float originalSize=0;
								if (Files.exists(filekk)){
									if (Files.isSymbolicLink(filekk)) {
										originalPath = Files
												.readSymbolicLink(filekk);
									}
									else
										originalPath=filekk;
									System.out.println("File:"+originalPath.toAbsolutePath());
									if (Files.exists(originalPath)){
								    	originalSize=Files.size(originalPath);
								    	// SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss");
								    	// originalFileDate=sdf.format(Files.getLastModifiedTime(originalPath));
								    	originalFileDate=Files.getLastModifiedTime(originalPath).toString();
									}
									else{
										
										String s="";
								        Process p;
								        String realName="";
							        	// System.out.println("readlink -f "+bamFiles[kk]);
							            p = Runtime.getRuntime().exec("readlink -f "+bamFiles[kk]);
							            
							            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
							            while ((s = stdInput.readLine()) != null) {
							                    realName=s;
							            }
							            p.waitFor();
							            // System.out.println ("exit: " + p.exitValue());
							            p.destroy();
								        
								        // System.out.println(s);
								        if (Files.exists(FileSystems.getDefault().getPath(realName))){
								        	originalSize=Files.size(FileSystems.getDefault().getPath(realName));
								        	originalPath=FileSystems.getDefault().getPath(realName);
								        	originalFileDate=Files.getLastModifiedTime(FileSystems.getDefault().getPath(realName)).toString();
								        }
								        else{
								        	originalSize=0;
								        	originalFileDate="Original file is not exisited";
								        }
								        
									}
								}
								else{
									originalSize=0;
									originalFileDate ="Not exisited";
									originalPath =filekk;	
								}
								String[] cc = bamFiles[kk].split("/");
								if ((cc.length != 10) && (cc.length != 11)) {
									writerErr
											.write("Error: "
													+ holder.get(i)
															.getAbsolutePath()
													+ " the file path has no 10 or 11 parts!");
									writerErr.newLine();
									parsingErr = true;
									break;
								}
								String platform = cc[4];
								if (cc.length==11)
								  platform=cc[5];
								
								if (originalPath.toString().toUpperCase().contains("NEXTSEQ"))
									platform="NextSeq";
								if (originalPath.toString().toUpperCase().contains("MISEQ"))
									platform="MiSeq";
									
								String flowcell = cc[7];
								if (cc.length==11)
									flowcell=cc[8];
								String bam = cc[9];
								if (cc.length==11)
									bam=cc[10];
								String[] dd = bam.split("_");
								
								// For debug temp
								writerErr
								.write("In"
										+ holder.get(i)
												.getAbsolutePath()
										+ ", processing the BAM file name ("+bam+")...");
					         	writerErr.newLine();
								if ((dd.length != 3) && (dd.length!=9)) {
									if ((dd.length==4) || (dd.length==10)) {
										if (dd[1].length()==1  && dd[1].matches("\\d+")) {
											writerErr
											.write("Warning: in"
													+ holder.get(i)
															.getAbsolutePath()
													+ " the BAM file name ("+bam+") has four or ten parts and the second part is not a digit number for four or ten parts!");
										}
									
									    else {	
											writerErr
													.write("Error: in"
															+ holder.get(i)
																	.getAbsolutePath()
															+ " the BAM file name ("+bam+") has no three or nine parts and the second part is not a digit number for four or ten parts!");
											writerErr.newLine();
											parsingErr = true;
											break;
									}
									}
								}
								
								// Below is modified to handle sample names in the format "SBXXXXX_1"
								// String sampleID = dd[0];
								String sampleID="";
								String[] partsOriginalPath=originalPath.toString().split("/");
								String[] partOriginalFileName=partsOriginalPath[partsOriginalPath.length-1].split("_");
								String index = dd[1];
								String lllane=dd[2];
								if (partOriginalFileName[1].length()==1 && partOriginalFileName[1].matches("\\d+")) {
									sampleID=partOriginalFileName[0]+"_"+partOriginalFileName[1];
								    index=dd[2];	
								    lllane=dd[3];
								}
							    else
							    	sampleID=partOriginalFileName[0];
								
								
								if (lllane.startsWith("L00")) {
									String lane = lllane.substring(3, 4);
									LaneSampleNode oneLane=null;
									if (originalDownsampledPath==null)
										oneLane = new LaneSampleNode(
												flowcell, lane, sampleID, group,
												index, platform,
												originalPath.toString(),"null",originalSize,downsampledSize,originalFileDate,downsampledFileDate,holder.get(i).getAbsolutePath(), downsampleRatio,0,null);
										
									
									else
									    oneLane = new LaneSampleNode(
											flowcell, lane, sampleID, group,
											index, platform,
											originalPath.toString(),originalDownsampledPath.toString(),originalSize,downsampledSize,originalFileDate,downsampledFileDate,holder.get(i).getAbsolutePath(),downsampleRatio, 0,null);
									LaneSampleNode storedLane = ll.head;
									boolean errFound = false;

									while (storedLane != null) {
										if (oneLane.compareTo(storedLane) == true) {
											writerErr
													.write("Error: "
															+ holder.get(i)
																	.getAbsolutePath()
															+ " has duplicated lane-level BAM file!");
											writerErr.newLine();
											errFound = true;
											break;
										}
										storedLane = storedLane.next;

									}
									if (errFound) {
										parsingErr = true;
										break;
									} else {
										ll.add(oneLane);
									}
								} else {
									writerErr
											.write("Error: "
													+ holder.get(i)
															.getAbsolutePath()
													+ " has format error in the lane-level BAM file name: (not L00X)!");
									writerErr.newLine();
									parsingErr = true;
									break;
								}
							}// for ends
						}
						if (parsingErr) {
							break;
						} else {
							if (lineDownsampled.startsWith("Downsampled"))
							   line = readerLog.readLine(); // Read the following
															// line after
															// "donwsampled"
						
						    else
							   line=lineDownsampled;
							if (line == null) {
								writerErr.write("Error: "
										+ holder.get(i).getAbsolutePath()
										+ ": input is not followed by output!");
								writerErr.newLine();
								parsingErr = true;
								break;

							} else {
								if (line.startsWith("Output:")) {
									readerLog.readLine();
									line = readerLog.readLine();
									if (line.contains("Cleaning")
											|| line.contains("cleaning")) {
										countCleaning++;
										if (countCleaning > 2) {
											writerErr
													.write("Error: "
															+ holder.get(i)
																	.getAbsolutePath()
															+ " has two merged BAM file!");
											writerErr.newLine();
											parsingErr = true;
											break;
										} else
											analysisIDs.add(ll);
									} else {
										if (!line.contains("Skipped")
												&& !line.contains("skipped")) {
											writerErr
													.write("Error: "
															+ holder.get(i)
																	.getAbsolutePath()
															+ ": Output is not followed by cleaning up or skipped!");
											writerErr.newLine();
											parsingErr = true;
											break;
										}
									}
								} else {
									writerErr
											.write("Error: "
													+ holder.get(i)
															.getAbsolutePath()
													+ " downsampled is not followed by output:!");
									writerErr.newLine();
									parsingErr = true;
									break;
								}
							}
							if (parsingErr)
								break;
						}
					}// Input block
					else {
						writerErr.write("Error: "
								+ holder.get(i).getAbsolutePath()
								+ " Building is not followed by Input:!");
						writerErr.newLine();
						parsingErr = true;
						break;
					}
					
				} // building block
			}// while parsing
			// System.out.println("next for");
			// if (i==0) break;
			readerLog.close();
			// readerLog2.close();
		}//for block
		writerErr.flush();
		writerErr.close();
		writerOut.write("AnalysisID"+"\t"+"Flowcell"+"\t"+"Group"+"\t"
				+"Index"+"\t"+"Platform"+"\t"+"Lane"+"\t"+"sample ID"+"\t"
				+"Input BAM file"+"\t"+"Modified Date"+"\t"+"Size"+"\t"
				+"Downsampled BAM file"+"\t"+"Downsample Ratio"
				+"\t"+"Log File Name");
		writerOut.newLine();
		for(int i=0;i<analysisIDs.size();i++){
			System.out.println("Writing:"+i+"/"+analysisIDs.size());
	    	String analysisID=analysisIDs.get(i).analysisID;
	    	LaneSampleNode sampleNodeofAnalysisID=analysisIDs.get(i).head;
	    	
	    	do{
	    		writerOut.write(analysisID+"\t"+sampleNodeofAnalysisID.flowcell+"\t"+sampleNodeofAnalysisID.group+"\t"
	    				+sampleNodeofAnalysisID.index+"\t"+sampleNodeofAnalysisID.instrumentID+"\t"+sampleNodeofAnalysisID.lane+"\t"+sampleNodeofAnalysisID.sampleID+"\t"
	    				+sampleNodeofAnalysisID.filePath+"\t"+sampleNodeofAnalysisID.fileModifiedDate+"\t"+sampleNodeofAnalysisID.fileSize+"\t"
	    				+sampleNodeofAnalysisID.downsampledFilePath+"\t"+sampleNodeofAnalysisID.downsampleRatio
	    				+"\t"+sampleNodeofAnalysisID.logFileName);
	    		LaneSampleNode otherNodes=sampleNodeofAnalysisID.next;
	    		while(otherNodes!=null){
	    			if (sampleNodeofAnalysisID.compareTo(otherNodes)==true){
	    				writerOut.write("\t"+"Warning:same to other BAM");
	    				break;
	    			}
	    			otherNodes=otherNodes.next;
	    		}
	    		writerOut.newLine();
	    		sampleNodeofAnalysisID=sampleNodeofAnalysisID.next;
	    	}
	    	while(sampleNodeofAnalysisID!=null);
	    }
		writerOut.flush();
		writerOut.close();
		System.out.println("Done!");
	}
}
