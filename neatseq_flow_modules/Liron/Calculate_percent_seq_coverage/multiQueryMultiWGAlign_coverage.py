#!/usr/bin/python
import subprocess
import csv
import os
from pprint import pprint
import sys
import argparse
import shlex


csv.field_size_limit(sys.maxsize)

def CommandLineProcessor():
	global samplesIndexFile
	global bpCoverageTreshold
	global outFile
	global targetSeqFile
	

	parser = argparse.ArgumentParser(description='cal_bam_percent_seq_coverage.py\n A pipeline to calculate the percent cover (~alignment) of each of a set of reference sequences within each sample\'s BAM file.')
	parser.add_argument('-s','--sindex', help='A file where each line defines one sample\'s name and its corresponding BAM file, tab-delimited.', required=True)
	parser.add_argument('-r','--reference_seqs', help='A multi-fasta file containing the refrences sequences whose sequence coverage we want to calculate.', required=True)
	parser.add_argument('-o','--output_file', help='A file listing the results (coverage of each sequence within each sample, in percents), rows are reference sequences and columns are samples, tab-delimited.', required=True)
	parser.add_argument('-t','--treshold', help='The minimum number of times a given sequence\'s basepair (position) need to be hit by the sample\'s reads, in order to be considered "a covered position".', required=True, type=int)
	
	#parser.add_argument('-t','--threads', help='Number of threads to be used for the bowtie2-align command.', required=True, type=int)


	
	args = parser.parse_args()
	samplesIndexFile = args.sindex
	bpCoverageTreshold = args.treshold
	outFile = args.output_file
	targetSeqFile = args.reference_seqs
	
	if not os.path.exists(targetSeqFile):
		print("Error: reference fasta file not found!")
		sys.exit(1)
	if not os.path.exists(samplesIndexFile):
		print("Error: samples index file not found!")
		sys.exit(1)
	# if os.path.exists(outFile):
		# print("Error: output file already exists!!")
		# sys.exit(1)


efh=open('err','wb')
sofh=open('out','wb')

def readSamplesIndexFile(fn):
	sampleIndexFile={}
	fh=open(fn, 'rb')
	csvreader=csv.reader(fh,delimiter="\t")
	for line in csvreader:
		if len(line)!=2:
			print("Error: all lines in samples index file must be tab-delimited with two items: sample_name, sample_file!")
			sys.exit(1)
		if not os.path.exists(line[1]):
			print("Error: the following BAM file does not exist - \"%s\"" % line[1])
			sys.exit(1)
		if line[0] in sampleIndexFile:
			print("Error: the following sample name appear twice: \"%s\"" % line[0])
		sampleIndexFile[line[0]]=line[1]
	fh.close()
	
	return sampleIndexFile
def countSeqLenInFasta(fastaFilename):
	seqlen=0
	seqname="**UNDEFINED**"
	fh=open(fastaFilename, 'rb')
	seqLenDict={}
	for line in fh:
		line=line.rstrip()
		if line[0]=='>':
			if seqname!="**UNDEFINED**":
				seqLenDict[seqname]=seqlen
			seqlen=0
			seqname=line[1:]
		elif seqname!="**UNDEFINED**":
			seqlen+=len(line.rstrip().replace('N',''))
	
	seqLenDict[seqname]=seqlen
	fh.close()
	
	return seqLenDict

	
	
def LaunchSubProcess(ProcParam, stdlogFH, errlogFH):
	print("Launching command: \"%s\"") % ProcParam
	proc=subprocess.Popen(ProcParam, shell=True,stdout=stdlogFH, stderr=errlogFH, universal_newlines=True)
	proc.wait()
	rc = proc.returncode
	if (rc != 0):
		print "Error: the program just executed has returned a non-zero exit code!"
		sys.exit()
#	return proc
	
def readCoveragePerBaseFromPILEUP(mpileupFile):
	print "Reading coverage info from the pile-up file \"%s\".." %mpileupFile
	mySequenceCoverage={}
	fh=open(mpileupFile, 'rb')
	csvreader=csv.reader(fh, delimiter="\t", quoting = csv.QUOTE_NONE, quotechar = "")
	for line in csvreader:
		seqName=line[0]
		baseCoverage=int(line[3])
		if not seqName in mySequenceCoverage:
			mySequenceCoverage[seqName]=[]
		mySequenceCoverage[seqName].append(baseCoverage)
	
	fh.close()

	return mySequenceCoverage

def countCoverageBasesCoveragWithMinTreshold(mySequencesCoverage, minBaseCoverage):
	CoverageBasesCoveragWithMinTresholdCountDict={}
	for seqName in mySequencesCoverage:
		CoverageBasesCoveragWithMinTresholdCountDict[seqName]=0
		for baseCoverage in mySequencesCoverage[seqName]:
			if (baseCoverage > minBaseCoverage):
				CoverageBasesCoveragWithMinTresholdCountDict[seqName]+=1
	return CoverageBasesCoveragWithMinTresholdCountDict


def calculateRelativeCoverageOfSequences(seqCoveredBasesCount,seqLengths):
	coverageAsPercent={}
	for seqname in seqLengths:
		if not seqname in seqCoveredBasesCount:
			coverageAsPercent[seqname]=0
		else:
			coverageAsPercent[seqname]=float(seqCoveredBasesCount[seqname])/float(seqLengths[seqname])*100.0
	return coverageAsPercent

def createPileupFile(bam_fn):
	mpileup=bam_fn+'.mpileup'
	if not os.path.exists(mpileup):
		print('Creating an mpileup file for the bam file \"%s\"..' % bam_fn)
		p=LaunchSubProcess('samtools mpileup --count-orphans -s ' + bam_fn + ' > ' + mpileup, sofh, efh)
#		p.wait()
	else:
		print('mpileup file for the bam file \"%s\" already exists!' % bam_fn)
		
def calSeqPercentCoverForBAMDict(bam_fn, seqLengths):
	global bpCoverageTreshold
	my_pileup=bam_fn+'.mpileup'
	createPileupFile(bam_fn)
	mySequencesCoverage=readCoveragePerBaseFromPILEUP(my_pileup)
	mySequencesBPCoveraedAboveTreshold=countCoverageBasesCoveragWithMinTreshold(mySequencesCoverage, bpCoverageTreshold)
	percentCoverPerSeq=calculateRelativeCoverageOfSequences(mySequencesBPCoveraedAboveTreshold,seqLengths)
	return percentCoverPerSeq

def writeResultsToFile(SequncePercentCoverInSamplesDictParam, outFile):
	samplesOrder=sorted(SequncePercentCoverInSamplesDictParam.keys())
	seqsOrder=sorted(SequncePercentCoverInSamplesDictParam[samplesOrder[0]].keys())
	fh=open(outFile, 'wb')
	csvwriter=csv.writer(fh, delimiter="\t")
	csvwriter.writerow(['Sequence/Sample'] + samplesOrder)

	for seq in seqsOrder:
		outline=[seq] + [SequncePercentCoverInSamplesDictParam[sample][seq] for sample in samplesOrder]
		csvwriter.writerow(outline)
	fh.close()

CommandLineProcessor()



seqLengths=countSeqLenInFasta(targetSeqFile)
SequncePercentCoverInSamplesDict={}
samplesIndex=readSamplesIndexFile(samplesIndexFile)

for sample in samplesIndex:
	print("Now processing the id \"%s\".." % sample)
	SequncePercentCoverInSamplesDict[sample]=calSeqPercentCoverForBAMDict(samplesIndex[sample], seqLengths)

efh.close()
sofh.close()
writeResultsToFile(SequncePercentCoverInSamplesDict, outFile)
#pprint(SequncePercentCoverInSamplesDict)
