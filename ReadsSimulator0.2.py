#!/bin/python

import sys, re, getopt, random, os, math
from Bio import SeqIO
#ys.path.append('./')
#import argument_process

# Parameter Initiation
#refGenome = sys.argv[1]
#readsOutFilePrefix = sys.argv[2]
#snpOutFile = sys.argv[3]
coverage = 50
errorRate = 0.0
readLen = 300
snpPercentage = 0.01
pairedEndReadLen = 100
snpFraction = 0.5
scoreRange=[0,40]
bufferRegion = 50
read1Len = pairedEndReadLen
read2Len = pairedEndReadLen
singleEndFlag = False
#seed = 10
#random.seed(seed)


# Classes
class Read:
    def __init__(self,readIndex,subSequence,idPrefix):
        self.startPosition = readIndex
        self.id = '@'+idPrefix+'_SEQID' + str(readIndex)
	#self.score = self.scoreGenerator(subSequence)
	self.scoreGenerator(subSequence)
        self.seq = subSequence
        self.seqLen = len(subSequence)
        if(errorRate > 0):
            self.addError()
        if not singleEndFlag:
            self.pairedEndConverter()        

    def addError(self):
        candAlleles = ['A','C','G','T']
        for i in range(self.seqLen):
            p = random.uniform(0,1)
            if(p <= errorRate):
                candAlleles.remove(self.seq[i])
                tempSeq = list(self.seq)
                tempSeq[i] = random.sample(candAlleles,1)[0]
                self.seq = ''.join(tempSeq)
                candAlleles = ['A','C','G','T'] 

    def scoreGenerator(self,seq):
        #rawScore = [scoreRange[1]]*len(seq)
        if singleEndFlag:
            rawScore = [int(scoreRange[1] - random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(self.seqLen)]
            self.score = ''.join([chr(x + 33) for x in rawScore])
        else:
            rawScore1 = [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(read1Len)]
            rawScore2 =  [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(read2Len)]
            self.read1Score = ''.join([chr(x+33) for x in rawScore1])
            self.read2Score = ''.join([chr(x+33) for x in rawScore2])

    def addSNP(self, snp):
        for index in range(len(snp.snpBufferIndex)):
            p = random.uniform(0,1)
            if(snp.snpBufferIndex[index] >= self.startPosition and snp.snpBufferIndex[index] < self.startPosition + readLen and p <= snp.snpFraction[index]):
                tempSeq = list(self.seq)
                tempSeq[snp.snpBufferIndex[index] - self.startPosition] = snp.snpAlleles[index] 
                self.seq = ''.join(tempSeq)
                snp.addedCount[index] = snp.addedCount[index] + 1
              
    def addIndel(self):
        pass
    
    def pairedEndConverter(self):
        self.read1 = self.seq[0:read1Len]
        baseComplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        tmpRead2 = self.seq[self.seqLen:self.seqLen-read2Len-1:-1]
        self.read2 = ''.join([baseComplement[x] for x in tmpRead2])
        #self.read1Score = self.score[0:read1Len]
        #self.read2Score = self.score[self.seqLen:self.seqLen-read2Len-1:-1]
        
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
    def __init__(self,snpPercentage,sequence,absoluteStart,strand,recordID):
        #random.seed(10) #delete later
        self.snpCount = int(math.floor(len(sequence) * snpPercentage))
        self.snpIndex = sorted(random.sample(range(len(sequence)),self.snpCount))
        self.snpBufferIndex = [x + bufferRegion for x in self.snpIndex]
        self.snpAlleles = self.findAllele(self.snpIndex,sequence)
        self.snpFraction = [snpFraction] * self.snpCount    
        self.addedCount = [0] * len(self.snpIndex)
        self.sequence = sequence
        self.recordID = recordID
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
            #random.seed(10) #delete later
            snpAlleles = snpAlleles + random.sample(candAlleles,1)
            candAlleles = ['A','C','G','T']
        return snpAlleles
    
    def writeToFile(self):
        with open(snpOutFile,'a') as f:
            for i in range(len(self.addedCount)):
                f.write('\t'.join([self.recordID,str(self.absolutePosition[i]),self.sequence[self.snpIndex[i]],self.snpAlleles[i],str(self.addedCount[i]),str(self.snpFraction[i])])+'\n')
        f.close()

    def printer(self):
        print self.snpIndex
	print self.snpAlleles
        print self.snpFraction	


class ReadsFileHandle:
    def __init__(self,filePrefix):
        if singleEndFlag:
            self.fileName = [filePrefix+'.fastq']
        else:
            self.fileName = [filePrefix+'1.fastq',filePrefix+'2.fastq']

    def fileOpen(self):
        if singleEndFlag:
            self.readsF = open(self.fileName[0],'a')
        else:
            self.reads1F = open(self.fileName[0],'a')
            self.reads2F = open(self.fileName[1],'a')
            
    def fileClose(self):
        if singleEndFlag:
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
def generateReads(seqBiopython,readsFileH):
    sequence = seqBiopython.seq
    idPrefix = seqBiopython.id
    geneStart = int(idPrefix.split('_')[3])
    strand = idPrefix.split('_')[1]
    snp = SNP(snpPercentage,sequence[bufferRegion:-bufferRegion],geneStart,strand,idPrefix)
    for i in range(len(sequence)-readLen):
        p = random.uniform(0,1)
	subSeq = sequence[i:i+readLen]
        if singleEndFlag:
            if(p <= coverage/float(readLen)):
	        read = Read(i,subSeq,idPrefix)
	        read.addSNP(snp)
                read.writeToFile(readsFileH)
        else:
            if(p <= coverage/float(read1Len+read2Len)):   
	        read = Read(i,subSeq,idPrefix)
	        read.addSNP(snp)
                read.pairedEndConverter()
                read.writeToFile(readsFileH)
    snp.writeToFile()

def argumentParser(argv):
    try:
        opts, args = getopt.getopt(argv,'hG:S:F:12',['help','coverage=','readsOutFile=','errorRate=','readLength=','snpFraction=','scoreRange=','bufferRegion=','snpPercentage=','snpOutFile=','refGenome=','readsOutFile='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt,arg in opts:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ('-G','--refGenome'):
            global refGenome
            refGenome = arg
        elif opt in ('-S','--snpOutFile'):
            global snpOutFile
            snpOutFile = arg
        elif opt in ('-F','--readsOutFile'):
            global readsOutFilePrefix
            readsOutFilePrefix = arg
        elif opt == '--errorRate':
            global errorRate
            errorRate = float(arg)
        elif opt == '--readLength':
            global readLength
            readLength = int(arg)
        elif opt == '--snpFraction':
            global snpFraction
            snpFraction = float(arg)
        elif opt == '--scoreRange':
            global scoreRange
            scoreRange = arg
        elif opt == '--bufferRegion':
            global bufferRegion
            bufferRegion = int(arg)
        elif opt == '--snpPercentage':
            global snpPercentage
            snpPercentage = float(arg)
        elif opt == '--coverage':
            global coverage
            coverage = int(arg)
        elif opt == '-1':
            global read1Len
            read1Len = int(arg)
        elif opt == '-2':
            global read2Len
            read2Len = int(arg)
        else:
            assert False, 'unhandeld option'
    if read1Len == 0 and read2Len == 0:
        singleEndFlag = True
    else:
        singleEndFlag = False

def usage():
    pass


def main():
    argumentParser(sys.argv[1:])
    readsFileH = ReadsFileHandle(readsOutFilePrefix)
    readsFileH.fileDelete()
    readsFileH.fileOpen()

    if os.path.isfile(snpOutFile):
        print 'Removing ' + snpOutFile + '!'
        os.system('rm ' + snpOutFile)

    for seqRecord in SeqIO.parse(refGenome,'fasta'):
        generateReads(seqRecord,readsFileH)

    readsFileH.fileClose()

if __name__ == "__main__":
    main()
