#!/bin/python

import sys
import re
from Bio import SeqIO
import random
import os

# Parameter Initiation
refGenome = sys.argv[1]
readsOutFilePrefix = sys.argv[2]
snpOutFile = sys.argv[3]
coverage = 50
errorRate = 0.0001
readLen = 300
snpCount = 20
pairedEndReadLen = 100
snpFraction = 0.5
scoreRange=[0,40]
bufferRegion = 50
read1Len = pairedEndReadLen
read2Len = pairedEndReadLen

# Classes
class Read:
    def __init__(self,readIndex,subSequence,idPrefix):
        self.startPosition = readIndex
        self.id = '@'+idPrefix+'_SEQID' + str(readIndex)
	self.score = self.scoreGenerator(subSequence)
	self.seq = subSequence
        self.seqLen = len(subSequence)

    def addError(self):
        candAlleles = ['A','C','G','T']
        for i in range(self.seqLen):
            p = random.uniform(0,1)
            if(p <= errorRate):
                candAlleles.remove(sequence[i])
                self.seq[i] = random.sample(candAlleles,1)
                candAlleles = ['A','C','G','T'] 

    def scoreGenerator(self,seq):
        rawScore = [scoreRange[1]]*len(seq)
        sangerScore = [chr(x+33) for x in rawScore]
        return "".join(sangerScore)

    def addSNP(self, snp):
        for index in range(len(snp.snpBufferIndex)):
            p = random.uniform(0,1)
            if(snp.snpBufferIndex[index] >= self.startPosition and snp.snpBufferIndex[index] < self.startPosition + readLen and p <= snp.snpFraction[index]):
                tempSeq = list(self.seq)
                tempSeq[snp.snpBufferIndex[index] - self.startPosition] = snp.snpAlleles[index] 
                self.seq = "".join(tempSeq)

    def addIndel(self):
        pass
    
    def pairedEndConverter(self):
        self.read1 = self.seq[0:read1Len]
        baseComplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        tmpRead2 = self.seq[self.seqLen:self.seqLen-read2Len-1:-1]
        self.read2 = ''.join([baseComplement[x] for x in tmpRead2])
        self.read1Score = self.score[0:read1Len]
        self.read2Score = self.score[self.seqLen:self.seqLen-read2Len-1:-1]
        
    def writeToFile(self,readsFileHandle):
        if read1Len == 0 and read2Len == 0:
            readsFileHandle.readsF.write(self.id + '\n')
            readsFileHandle.readsF.write(str(self.seq) + '\n')
            readsFileHandle.readsF.write('+' + '\n')
            readsFileHandle.readsF.write(self.score + '\n')
        elif read1Len > 0 and read2Len > 0:
            readsFileHandle.reads1F.write(self.id + '/1\n')
            readsFileHandle.reads1F.write(str(self.read1) + '\n')
            readsFileHandle.reads1F.write('+' + '\n')
            readsFileHandle.reads1F.write(self.read1Score + '\n')
            
            readsFileHandle.reads2F.write(self.id + '/2\n')
            readsFileHandle.reads2F.write(str(self.read2) + '\n')
            readsFileHandle.reads2F.write('+' + '\n')
            readsFileHandle.reads2F.write(self.read2Score + '\n')
        else:
            print 'Error: Check the read length!'
    def printer(self):
        print self.id
        print self.seq
        print self.score


class SNP:
    def __init__(self,snpCount,sequence,absoluteStart,strand):
        random.seed(10) #delete later
        self.snpIndex = sorted(random.sample(range(len(sequence)),snpCount))
        self.snpBufferIndex = [x + bufferRegion for x in self.snpIndex]
        self.snpAlleles = self.findAllele(self.snpIndex,sequence)
        self.snpFraction = [snpFraction] * snpCount    
        self.sequence = sequence
        self.strand = strand
        if strand == '+': 
            self.absolutePosition = [x + absoluteStart for x in self.snpIndex]
        elif strand == '_':
            self.absolutePosition = [absoluteStart - x + len(self.sequence)-1 for x in self.snpIndex]
        else:
            print 'Error: Check the strand!'

    def findAllele(self,snpIndex,sequence):
        snpAlleles = []
        candAlleles = ['A','C','G','T']
	for x in self.snpIndex:
            candAlleles.remove(sequence[x])
            random.seed(10) #delete later
            snpAlleles = snpAlleles + random.sample(candAlleles,1)
            candAlleles = ['A','C','G','T']
        return snpAlleles
    
    def writeToFile(self):
        with open(snpOutFile,'a') as f:
            for i in range(snpCount):
                f.write('\t'.join([seqRecord.id,str(self.absolutePosition[i]),self.sequence[self.snpIndex[i]],self.snpAlleles[i],str(self.snpFraction[i])])+'\n')
        f.close()

    def printer(self):
        print self.snpIndex
	print self.snpAlleles
        print self.snpFraction	


class ReadsFileHandle:
    def __init__(self,filePrefix):
        if read1Len == 0 and read2Len == 0:
            self.fileName = [filePrefix+'.fastq']
        else:
            self.fileName = [filePrefix+'1.fastq',filePrefix+'2.fastq']

    def fileOpen(self):
        if read1Len == 0 and read2Len ==0:
            self.readsF = open(self.fileName[0],'a')
        else:
            self.reads1F = open(self.fileName[0],'a')
            self.reads2F = open(self.fileName[1],'a')
            
    def fileClose(self):
        if read1Len == 0 and read2Len ==0:
            self.readsF.close()
        else:
            self.reads1F.close()
            self.reads2F.close()
    
    def fileDelete(self):
        for file in self.fileName:
            if os.path.isfile(file):
                print 'Removing ' + str(file) + '!'
                os.system('rm ' + file)

# Function for generating reads
def generateReads(seqBiopython):
    sequence = seqBiopython.seq
    idPrefix = seqBiopython.id
    geneStart = int(idPrefix.split('_')[3])
    strand = idPrefix.split('_')[1]
    snp = SNP(snpCount,sequence[bufferRegion:-bufferRegion],geneStart,strand)
    snp.writeToFile()
    for i in range(len(sequence)-readLen):
        p = random.uniform(0,1)
	subSeq = sequence[i:i+readLen]
        if read1Len == 0 and read2Len == 0:
            if(p <= coverage/float(readLen)):
	        read = Read(i,subSeq,idPrefix)
                #read.AddError()
	        read.addSNP(snp)
                read.writeToFile(readsFileH)
        elif read1Len > 0 and read2Len >0:
            if(p <= coverage/float(read1Len+read2Len)):   
	        read = Read(i,subSeq,idPrefix)
	        read.addSNP(snp)
                read.pairedEndConverter()
                read.writeToFile(readsFileH)
        else:
            exit('Error: Check the read length!')

readsFileH = ReadsFileHandle(readsOutFilePrefix)
readsFileH.fileDelete()
readsFileH.fileOpen()

if os.path.isfile(snpOutFile):
    print 'Removing ' + snpOutFile + '!'
    os.system('rm ' + snpOutFile)

for seqRecord in SeqIO.parse(sys.argv[1],'fasta'):
    generateReads(seqRecord)

readsFileH.fileClose()


