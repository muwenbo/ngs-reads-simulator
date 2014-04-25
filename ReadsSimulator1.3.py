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
readLen = 400
snpPercentage = 0.01
pairedEndReadLen = 100
snpFraction = 0.5
scoreRange=[0,40]
bufferRegion = 100
read1Len = pairedEndReadLen
read2Len = pairedEndReadLen
singleEndFlag = False
#seed = 10
#random.seed(seed)


# Classes
class Read:
    def __init__(self,readIndex,subSequence,idPrefix,error):
        self.startPosition = readIndex
        self.id = '@'+idPrefix+'_SEQID' + str(readIndex)
	#self.score = self.scoreGenerator(subSequence)
	self.scoreGenerator(subSequence)
        self.seq = subSequence
        self.seqLen = len(subSequence)
        self.suffix = ''
        if(errorRate > 0):
            self.addError(error)
        if not singleEndFlag:
            self.pairedEndConverter()        

    def addError(self,error):
        candAlleles = ['A','C','G','T']
        for i in range(self.seqLen):
            p = random.uniform(0,1)
            if(p <= errorRate):
                candAlleles.remove(self.seq[i])
                tempSeq = list(self.seq)
                rawAllele = tempSeq[i]
                tempSeq[i] = random.sample(candAlleles,1)[0]
                self.seq = ''.join(tempSeq)
                candAlleles = ['A','C','G','T'] 
                error.addItem(self.id,self.startPosition,i,rawAllele,tempSeq[i])
                if singleEndFlag:
                    self.score[i] = chr(random.sample(range(0,20),1)[0]+33)
                else:
                    if i < read1Len:
                        self.read1Score[i] = chr(random.sample(range(0,20),1)[0]+33)
                    elif i >= self.seqLen - read2Len:
                        self.read2Score[i-(self.seqLen-read2Len)] = chr(random.sample(range(0,20),1)[0]+33)

    def scoreGenerator(self,seq):
        #rawScore = [scoreRange[1]]*len(seq)
        if singleEndFlag:
            rawScore = [int(scoreRange[1] - random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(self.seqLen)]
            self.score = [chr(x + 33) for x in rawScore]
        else:
            rawScore1 = [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(read1Len)]
            rawScore2 =  [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in range(read2Len)]
            self.read1Score = [chr(x+33) for x in rawScore1]
            self.read2Score = [chr(x+33) for x in rawScore2]

    def addSNP(self, snp):
        for index in range(len(snp.snpBufferIndex)):
            p = random.uniform(0,1)
            if(singleEndFlag):
                if(snp.snpBufferIndex[index] >= self.startPosition and snp.snpBufferIndex[index] < self.startPosition + self.seqLen):
                    snp.totalReadsCount[index] = snp.totalReadsCount[index] + 1
                    if(p <= snp.snpFraction[index]+0.01):
                        tempSeq = list(self.seq)
                        tempSeq[snp.snpBufferIndex[index] - self.startPosition] = snp.snpAlleles[index] 
                        self.seq = ''.join(tempSeq)
                        snp.addedCount[index] = snp.addedCount[index] + 1
                        self.suffix = self.suffix + '_' + str(snp.absolutePosition[index])
            else:
                 if((snp.snpBufferIndex[index] >= self.startPosition and snp.snpBufferIndex[index] < self.startPosition + read1Len) or (snp.snpBufferIndex[index] >= self.startPosition + self.seqLen - read2Len and snp.snpBufferIndex[index] < self.startPosition + self.seqLen)):
                    snp.totalReadsCount[index] = snp.totalReadsCount[index] + 1
                    if(p <= snp.snpFraction[index]+0.01):
                        tempSeq = list(self.seq)
                        tempSeq[snp.snpBufferIndex[index] - self.startPosition] = snp.snpAlleles[index] 
                        self.seq = ''.join(tempSeq)
                        snp.addedCount[index] = snp.addedCount[index] + 1
                        self.suffix = self.suffix + '_' + str(snp.absolutePosition[index])
                              

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
            readsFileHandle.readsF.write(self.id + self.suffix + '\n')
            readsFileHandle.readsF.write(str(self.seq) + '\n')
            readsFileHandle.readsF.write('+' + '\n')
            readsFileHandle.readsF.write(''.join(self.score) + '\n')
        elif read1Len > 0 and read2Len > 0:
            readsFileHandle.reads1F.write(self.id + self.suffix + '/1\n')
            readsFileHandle.reads1F.write(str(self.read1) + '\n')
            readsFileHandle.reads1F.write('+' + '\n')
            readsFileHandle.reads1F.write(''.join(self.read1Score) + '\n')
            
            readsFileHandle.reads2F.write(self.id + self.suffix + '/2\n')
            readsFileHandle.reads2F.write(str(self.read2) + '\n')
            readsFileHandle.reads2F.write('+' + '\n')
            readsFileHandle.reads2F.write(''.join(self.read2Score) + '\n')
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
        self.totalReadsCount = [0] * len(self.snpIndex)
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
    
    def writeToFile(self,fileH):
        for i in range(len(self.addedCount)):
            fileH.write('\t'.join([self.recordID,str(self.absolutePosition[i]),self.sequence[self.snpIndex[i]],self.snpAlleles[i],str(self.addedCount[i]),str(self.totalReadsCount[i])])+'\n')

    def printer(self):
        print self.snpIndex
	print self.snpAlleles
        print self.snpFraction	


class Error:
    def __init__(self):
        self.errorInfor = {}

    def addItem(self,id,parentIndex,localIndex,rawAllele,newAllele):
        globalIndex = str(parentIndex + localIndex)
        if globalIndex not in self.errorInfor:
            self.errorInfor[globalIndex] = [id,rawAllele,newAllele]
        else:
            self.errorInfor[globalIndex][2] = self.errorInfor[globalIndex][2] + '\\' + newAllele    

    def writeToFile(self,fileH):
        for key in self.errorInfor:
            fileH.write('\t'.join(self.errorInfor[key])+ '\t' + key + '\n')


class ReadsFileHandle:
    def __init__(self,filePrefix):
        if singleEndFlag:
            self.fileName = [filePrefix+'.fastq']
        else:
            self.fileName = [filePrefix+'1.fastq',filePrefix+'2.fastq']

    def fileOpen(self):
        self.errorFH = open(errorOutFile,'w')
        self.snpFH = open(snpOutFile,'w')
        if singleEndFlag:
            self.readsF = open(self.fileName[0],'a')
        else:
            self.reads1F = open(self.fileName[0],'a')
            self.reads2F = open(self.fileName[1],'a')
            
    def fileClose(self):
        self.errorFH.close()
        self.snpFH.close()
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
        if os.path.isfile(snpOutFile):
            print 'Removing ' + snpOutFile + '!'
            os.system('rm ' + snpOutFile)

        if os.path.isfile(errorOutFile):
            print 'Removing ' + errorOutFile + '!'
            os.system('rm ' + errorOutFile)



# Function for generating reads
def generateReads(seqBiopython,readsFileH):
    sequence = seqBiopython.seq
    idPrefix = seqBiopython.id
    geneStart = int(idPrefix.split('_')[3])
    strand = idPrefix.split('_')[1]
    snp = SNP(snpPercentage,sequence[bufferRegion:-bufferRegion],geneStart,strand,idPrefix)
    error = Error()
    seqLen =len(sequence)
#    for i in range(seqLen - min(readLen,read1Len+read2Len)):
    for i in range(seqLen - readLen):
        if singleEndFlag:
            multiple = coverage/float(readLen)
            while(multiple > 0):
                p = random.uniform(0,1)
                if(p <= coverage/float(readLen)):
	            subSeq = sequence[i:i+readLen]
	            read = Read(i,subSeq,idPrefix,error)
	            read.addSNP(snp)
                    read.writeToFile(readsFileH)
                multiple = multiple - 1
        else:
            multiple = coverage/float(read1Len+read2Len)
            while(multiple > 0):
                p = random.uniform(0,1)
                if(p <= multiple):   
	            #dynamicGap = random.sample(range(min(readLen,seqLen - i)-read1Len-read2Len),1)[0]
                    dynamicGap = max(random.sample(range(200,300),1)[1],readLen-read1Len,read2Len)
                    newReadLen = read1Len + dynamicGap + read2Len
                    subSeq = sequence[i:i + newReadLen]
	            read = Read(i,subSeq,idPrefix,error)
	            read.addSNP(snp)
                    read.pairedEndConverter()
                    read.writeToFile(readsFileH)
                multiple = multiple-1
    snp.writeToFile(readsFileH.snpFH)
    error.writeToFile(readsFileH.errorFH)

def argumentParser(argv):
    try:
        opts, args = getopt.getopt(argv,'hG:S:F:E:12',['help','coverage=','readsOutFile=','errorRate=','readLength=','snpFraction=','scoreRange=','bufferRegion=','snpPercentage=','snpOutFile=','refGenome=','readsOutFile=','errorOutFile'])
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
        elif opt in ('-E','--errorOutFile'):
            global errorOutFile
            errorOutFile = arg
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
    print 'Start simulating.'
    for seqRecord in SeqIO.parse(refGenome,'fasta'):
        generateReads(seqRecord,readsFileH)

    readsFileH.fileClose()

if __name__ == "__main__":
    main()
