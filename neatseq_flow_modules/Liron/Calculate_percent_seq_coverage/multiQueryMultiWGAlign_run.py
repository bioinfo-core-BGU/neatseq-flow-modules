#!/usr/bin/python
from collections import deque
import subprocess;
import sys
import time
import argparse
import csv
import os

def CommandLineProcessor():
	global genomesListFiles
	global readsFilesListFile
	global threads_limit
	global commandLine
	global simulate
	global outDir
	global stdOutFile
	global bowtie2_align_processes_limit
	global bowtie2_params
	global readsAreUnpaired
	global maxAlignsPerRead

	parser = argparse.ArgumentParser(description='multiQueryMultiWGAlign_run.py\n A pipeline to allign multiple queries to multiple whole genomes using bowtie2.')
	parser.add_argument('-g','--genomes', help='A file listing the genomes to be aligned, each line containing genome_id and path to the query fasta file, tab-delimited.', required=True)
	parser.add_argument('-r','--reads', help='A file listing the read files to be aligned, each line containing sample_id and path to the query R1 fastq file and R2 fastq file, tab-delimited.', required=True)
	parser.add_argument('-t','--threads', help='Number of threads to be used for the bowtie2-align command.', required=True, type=int)
	parser.add_argument('-p','--processes', help='Maximum number of bowtie2 alignment processes to run at one.', required=True, type=int)
	parser.add_argument('-m','--max_hits', help='Maximum number of alignments to allow per read', required=False, type=int, default = 1)
	#parser.add_argument('-s','--simulate', help='Don\'t actually run the commands but just print them.', required=False,action="store_true")
	parser.add_argument('-o','--outdir', help='The name of the output directory for the analysis.', required=True)
	parser.add_argument('-u','--unpaired', help='Tell bowtie to treat the reads as singles (unpaired)', required=False, action="store_true")
	#parser.add_argument('-b','--bowtie2_params', help='Comman line parameters to pass to bowtie2,', required=False)


	
	args = parser.parse_args()
	readsAreUnpaired = args.unpaired
	bowtie2_align_processes_limit = args.processes
	genomesListFiles = args.genomes
	readsFilesListFile = args.reads
	outDir = args.outdir
	maxAlignsPerRead = args.max_hits
	#bowtie2_params = args.bowtie2_params
	
	threads_limit = args.threads

	#simulate = args.simulate
	
	if (threads_limit < 1):
		print "Error: threads must be > 1"
		sys.exit(1)

	if (bowtie2_align_processes_limit < 1):
		print "Error: processes must be > 1"
		sys.exit(1)


def CheckCompulsoryFile(Filename, Description):
	if not os.path.isfile(Filename):
		print "Compulsory " + Description + " file \"" + Filename + "\" not exists!"
		sys.exit(1)

def CheckIfMustFilesExist():
	global genomesListFiles
	global readsFilesListFile
	global outDir

	CheckCompulsoryFile(genomesListFiles, "list of genomes to align in fasta format")
	CheckCompulsoryFile(readsFilesListFile, "list of samples to align in fastq format")
	
	if os.path.isdir(outDir):
		print "Error: output directory already exists!"
		sys.exit(1)


def isFilenameLegal(fn):
	legalChars='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-'
	if len(fn)==0:
		return False
	for c in fn:
		if not c in legalChars:
			return False
	return True

def readSamplesTable(fn):
	samplesTable=[]
	fh=open(fn,'rb')
	csvreader=csv.reader(fh,delimiter="\t")
	for line in csvreader:
		samplesTable.append(line)
		if not isFilenameLegal(line[0]):
			print "Error: sample index must be a legal filename but it is not: \"" + line[0] + "\"!"
			sys.exit(1)
	fh.close()
	
	return samplesTable

def readGenomesTable(fn):
	genomesTable=[]
	fh=open(fn,'rb')
	csvreader=csv.reader(fh,delimiter="\t")
	for line in csvreader:
		genomesTable.append(line)
		if not isFilenameLegal(line[0]):
			print "Error: genome index must be a legal filename but it is not: \"" + line[0] + "\"!"
			sys.exit(1)	
	fh.close()
	
	return genomesTable

def calcSamplesDepthFromFastqFiles(samplesTable):
	
	samplesDepthDict = {}
	
	print("Calculating samples depth..")
	i=0
	for sample in samplesTable:
		i+=1
		print("Calculating depth of sample %s/%s" % (i,len(samplesTable)))
		fh=open(sample[1])
		lc=0
		for line in fh:
			lc+=1
		
		fh.close()
		samplesDepthDict[sample[0]] = lc
	return samplesDepthDict

def constructBowtie2BuildIndexCommand(genomeFile, genome_id, outDir):
	genomeFileBaseName = os.path.basename(genomeFile)
	#cmd = "/programs/bowtie2/bowtie2-2.2.9/bowtie2-build -f " + genomeFile + " " + outDir + "/genome_indices/" + genome_id
	cmd = "bowtie2-build -f " + genomeFile + " " + outDir + "/genome_indices/" + genome_id
	return cmd
	
def constructBowtie2AlignCommand(genome_id, sample_id, sample_R1_file, sample_R2_file, outDir, ncpu, readsAreUnpaired):
	bowtie2_align_outputfile = outDir + "/" + sample_id  + "/"  + sample_id + "_vs_" + genome_id
	genome_index_full_path = outDir + "/genome_indices/" + genome_id
	PairedReadFilenamesParam =  " -1 " +  sample_R1_file + " -2 " + sample_R2_file + " "
	unPairedReadFilenamesParam =  " -U " +  sample_R1_file + "," + sample_R2_file + " "
	if readsAreUnpaired:
		ReadFilenamesParam=unPairedReadFilenamesParam
	else:
		ReadFilenamesParam=PairedReadFilenamesParam
	
	#cmd = "bowtie2 " + ReadFilenamesParam + " -x " + genome_index_full_path + " -S " + bowtie2_align_outputfile + " -p " + str(ncpu) + " -k " + str(maxAlignsPerRead)  + " --no-una && samtools view -bS " + bowtie2_align_outputfile + " | samtools sort - -o " + bowtie2_align_outputfile + ".bam && rm " + bowtie2_align_outputfile
	cmd = "bowtie2 " + ReadFilenamesParam + " -x " + genome_index_full_path + " -S " + bowtie2_align_outputfile + " -p " + str(ncpu) + " -k " + str(maxAlignsPerRead)  + " --no-una;   samtools view -bS " + bowtie2_align_outputfile + " >  " + bowtie2_align_outputfile + ".sam" + " ; " + "samtools sort " + bowtie2_align_outputfile + ".sam" + " -o " + bowtie2_align_outputfile + ".bam  ; rm " + bowtie2_align_outputfile + " " + bowtie2_align_outputfile + ".sam"
	#samFile = bowtie2_align_outputfile + ".sam"
	#cmd1 = "bowtie2 " + ReadFilenamesParam + " -x " + genome_index_full_path + " -S " + bowtie2_align_outputfile + " -p " + str(ncpu) + " -k " + str(maxAlignsPerRead)  + " --no-una > " + samFile
	#cmd2 = "samtools view -bS " + samFile + " > " + bowtie2_align_outputfile
	#cmd3 = "samtools sort " + bowtie2_align_outputfile + " -o " + bowtie2_align_outputfile + ".bam"
	#cmd4 = "rm " + bowtie2_align_outputfile + " " + samFile
	
	#return [cmd1, cmd2, cmd3, cmd4]
	return cmd
	
def createBowtie2GenomeIndices(genomesTable, outDir):
	print("\nStage 1: creating genome indices..\n")
	genomesCount=len(genomesTable)
	i=0
	for genomeInfo in genomesTable:
		i+=1
		print "Creating index for genome %d/%d, genome_id: %s, genome FASTA file: \"%s\"" % (i,genomesCount,genomeInfo[0],genomeInfo[1])
		cmd = constructBowtie2BuildIndexCommand(genomeInfo[1], genomeInfo[0], outDir)
		print cmd
		errlogfile = open(outDir+'/genome_indices/error_log.txt','ab').close()
		stdlogfile = open(outDir+'/genome_indices/standard_log.txt','ab').close()
		errlogfile = open(outDir+'/genome_indices/error_log.txt','ab')
		stdlogfile = open(outDir+'/genome_indices/standard_log.txt','ab')
		subprocess.call(cmd, shell=True, stdout=stdlogfile, stderr=errlogfile, universal_newlines=True)
		errlogfile.close()
		stdlogfile.close()
	

def SubmitAlignmentProcessesGradually(samplesTable,genomesTable,outDir,nthreads,ncpu,readsAreUnpaired):
	subProcessesPool=[]
	subProcessesIDs=[]
	subProcessesSample_ids=[]
	activeSubProcCount=0
	JobsProcessedCount=0

	openSlots=nthreads

	errlogfile = open(outDir+'/error_log.txt','ab').close()
	stdlogfile = open(outDir+'/standard_log.txt','ab').close()
	errlogfile = open(outDir+'/error_log.txt','ab')
	stdlogfile = open(outDir+'/standard_log.txt','ab')


	print("\nStage 2: aligning samples to genomes..\n")
	
	for gindex in range(0,len(genomesTable)):
		for sindex in range(0,len(samplesTable)):
			#errlogfile = open(outDir+'/'+samplesTable[sindex][0]+'/error_log.txt','ab')
			#stdlogfile = open(outDir+'/'+samplesTable[sindex][0]+'/standard_log.txt','ab')
			print "Now aligning sample_id %s to genomoe_id %s" % (samplesTable[sindex][0], genomesTable[gindex][0])
			genomeIndexPath=outDir+"/genome_indices/"+genomesTable[gindex][0]
	
			cmd = constructBowtie2AlignCommand(genomesTable[gindex][0], samplesTable[sindex][0], samplesTable[sindex][1], samplesTable[sindex][2], outDir, ncpu, readsAreUnpaired)
			print(cmd)
			openSlots=0
			while (openSlots<1):
				activeSubProcCount=CountLiveSubProcesses(subProcessesPool)
				openSlots=nthreads - activeSubProcCount
				time.sleep(1)
			proc=LaunchSubProcess(cmd,stdlogfile, errlogfile)
			subProcessesPool.append(proc)
			subProcessesSample_ids.append(samplesTable[sindex][0])

	print "\nAll alignment jobs submitted!\n"
	activeSubProcCount=1
	while (activeSubProcCount>0):
				activeSubProcCount=CountLiveSubProcesses(subProcessesPool)
				time.sleep(1)

	errlogfile.close()
	stdlogfile.close()
			
	#for i in range(len(subProcessesSample_ids)):
	#		errlogfile=outDir+'/'+subProcessesSample_ids[i]+"/error.log"
	#		stdlogfile=outDir+'/'+subProcessesSample_ids[i]+"/stdout.log"
	#		errlogFH=open(errlogfile, 'ab')
	#		stdlogFH=open(stdlogfile, 'ab')
	#		out, err = proc.communicate()
	#		errlogFH.write(err)
	#		stdlogFH.write(out)
	#		errlogFH.close()
	#		stdlogFH.close()


def CountLiveSubProcesses(SubProcPool):
	n=0
	for proc in SubProcPool:
		if proc.poll()==None:
			n+=1
	return n

def LaunchSubProcess(ProcParam, stdlogFH, errlogFH):
	proc=subprocess.Popen(ProcParam, shell=True,stdout=stdlogFH, stderr=errlogFH, universal_newlines=True)
	return proc

def writeTSVFile(list2d,outFileame):
	fh=open(outFileame,'wb')
	csvwriter=csv.writer(fh,delimiter="\t",quoting = csv.QUOTE_NONE)
	for line in list2d:
		csvwriter.writerow(line)
	
	fh.close()

def writeTSVFileFromDict(mydict,outFileame):
	fh=open(outFileame,'wb')
	csvwriter=csv.writer(fh,delimiter="\t",quoting = csv.QUOTE_NONE)
	for k in mydict:
		csvwriter.writerow([k,mydict[k]])
	
	fh.close()

def WriteGenomesAndSamplesMetadata(samplesTable, genomesTable, outDir):
	writeTSVFile(samplesTable, outDir+'/samples_info.tsv')
	writeTSVFile(genomesTable, outDir+'/genomes_info.tsv')

def createSampleDirectories(samplesTable, outDir):
	print "Creating sampels sub-directories.."
	for sample in samplesTable:
		os.makedirs(outDir+"/"+sample[0])

CommandLineProcessor()
CheckIfMustFilesExist()
os.makedirs(outDir)
os.makedirs(outDir + "/genome_indices")
samplesTable=readSamplesTable(readsFilesListFile)
genomesTable=readGenomesTable(genomesListFiles)
createSampleDirectories(samplesTable, outDir)
samplesDepthDict=calcSamplesDepthFromFastqFiles(samplesTable)
writeTSVFileFromDict(samplesDepthDict,outDir+"/samples_depth.tsv")
WriteGenomesAndSamplesMetadata(samplesTable, genomesTable, outDir)
createBowtie2GenomeIndices(genomesTable, outDir)
SubmitAlignmentProcessesGradually(samplesTable,genomesTable,outDir,bowtie2_align_processes_limit,threads_limit,readsAreUnpaired)
