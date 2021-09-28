package nci.cgr.manifest;

import java.util.Date;

public class FastqFile {

	public FastqFile(String platform,String flowcell, String lane,  String index, String project,String sample,Date updateDate) {
		super();
		this.flowcell = flowcell;
		this.index = index;
		this.project = project;
		this.lane = lane;
		this.sample = sample;
		this.platform=platform;
		this.updateDate=updateDate;
	}
	
	String flowcell;
	String index;
	String project;
	String lane;
	String sample;
	String platform;
	Date updateDate;
	
	public String getFlowcell() {
		return flowcell;
	}
	public void setFlowcell(String flowcell) {
		this.flowcell = flowcell;
	}
	public String getIndex() {
		return index;
	}
	public void setIndex(String index) {
		this.index = index;
	}
	public String getProject() {
		return project;
	}
	public void setProject(String project) {
		this.project = project;
	}
	public String getLane() {
		return lane;
	}
	public void setLane(String lane) {
		this.lane = lane;
	}
	public String getSample() {
		return sample;
	}
	public void setSample(String sample) {
		this.sample = sample;
	}
	public String getPlatform() {
		return platform;
	}
	public void setPlatform(String platform) {
		this.platform = platform;
	}
	
	
}
