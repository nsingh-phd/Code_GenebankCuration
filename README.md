## This repository contains the `R` code and other required files to replicate the analysis performed in this manuscript.

1. Using the provided keyfiles, run Tassel4 pipeline to get hap files.
2. Run `1_SNPFiltering.R` file to filter and merge three hap files.
3. Run `2_AnalysisPipeline.R` file to perform the analysis.

### Set up Tassel4 and GBS pipeline in Eclipse
	* Requirements
		1. Eclipse - Version: Oxygen.3a Release (4.7.3a)
		2. Download `Genebank.java` to your computer
		3. Java 7
	* Follow these steps-
		1. Open Eclipse and set workspace
		2. Select File -> Import... -> General
		3. Select 'Projects from Folder or Archive'
		4. Click 'Directory' and select `Genebank.java` and click Open
		5. Select all and click Finish
	* You should have three projects imported in Eclipse
		* Genebank.java, PipelineTassel4, and tassel4
	* If you have errors (red x) in front of tassel4 or any other project then follow these steps
		1. Open Eclipse preferences
		2. Go to Java -> Installed JREs
		3. Select version 7 (also listed as 1.7)
		4. Apply

### Create runnable jar file
	* Open PipelineTassel4 -> src -> pipeline -> Tassel4Pipeline.java
	* Update variables - projectName, userName, keyFile etc.
	* Save and Run
	* Ignore the warnings if any
	* Right click on `Tassel4Pipeline.java` -> Export...
	* Select Java -> Runnable JAR File -> Next
	* In configuration select 'Tassel4Pipeline - PipelineTassel4'
	* Select Export destination
	* Select `Package required libraries...` under library handling
	* Click Finish and OK
	* Use the JAR file and Keyfile for SNP calling