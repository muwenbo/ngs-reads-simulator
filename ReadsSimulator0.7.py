#!/bin/python

import sys, re, getopt, random, os, math, array
from Bio import SeqIO

# Classes
class Read:
    def __init__(self,readIndex,sequence,idPrefix,error,scoreRange,errorRate,read1Len,read2Len):
        self.startPosition = readIndex
        self.id = '@'+idPrefix+'_SEQID' + str(readIndex)
        self.read1Len = read1Len
        self.read2Len = read2Len
        self.errorRate = errorRate
    #self.score = self.scoreGenerator(subSequence)
        dynamicGap = min(random.sample(range(200,300),1)[0],len(sequence)-readIndex-read1Len-read2Len)
        newReadLen = read1Len + dynamicGap + read2Len
        subSequence = sequence[readIndex:readIndex + newReadLen]                
        self.seq = array.array('c',subSequence).tolist()
        self.seqLen = len(subSequence)
        self.scoreGenerator(subSequence,scoreRange)
        self.suffix = ''
        if(errorRate > 0):
            self.addError(error)

    def addError(self,error):
        candAlleles = ['A','C','G','T']
        errorProb = [random.uniform(0,1) for x in range(self.seqLen)]
        errorIndex = [errorProb.index(x) for x in errorProb if x <= self.errorRate]    
        for i in errorIndex:
            rawAllele = self.seq[i]
            newAllele = random.sample([x for x in candAlleles if x != rawAllele],1)[0]
            self.seq[i] = newAllele
            error.addItem(self.id,self.startPosition,i,rawAllele,newAllele)
            self.score[i] = chr(random.sample(range(0,20),1)[0]+33)

    def scoreGenerator(self,seq,scoreRange):
        scoreR = range(int(math.ceil(self.seqLen/2.0)))
        rawScore1 = [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in scoreR]
        scoreR.reverse()
        if (self.seqLen/2)%1 != 0:
            scoreR = scoreR[:-1]
        rawScore2 = [int(scoreRange[1] -  random.gammavariate(1,0.5) - x*random.normalvariate(0.1,0.03)) for x in scoreR]
        self.score = [chr(x+33) for x in rawScore1 + rawScore2]
        #self.read1Score = [chr(x+33) for x in rawScore1]
        #self.read2Score = [chr(x+33) for x in rawScore2]

    def addSNP(self, snp):
        for index in range(len(snp.snpBufferIndex)):
            p = random.uniform(0,1)
            if((snp.snpBufferIndex[index] >= self.startPosition and snp.snpBufferIndex[index] < self.startPosition + self.read1Len) or (snp.snpBufferIndex[index] >= self.startPosition + self.seqLen - self.read2Len and snp.snpBufferIndex[index] < self.startPosition + self.seqLen)):
                snp.totalReadsCount[index] = snp.totalReadsCount[index] + 1
                if(p <= snp.snpFraction[index]+0.01):
                    self.seq[snp.snpBufferIndex[index] - self.startPosition] = snp.snpAlleles[index] 
                    snp.addedCount[index] = snp.addedCount[index] + 1
                    self.suffix = self.suffix + '_' + str(snp.absolutePosition[index])
                              

    def addDeletion(self,dels):
        for index in range(len(dels.delBufferIndex)):
            p = random.uniform(0,1)
            if((dels.delBufferIndex[index]+dels.delLen[index] >= self.startPosition and dels.delBufferIndex[index] < self.startPosition + self.read1Len) or (dels.delBufferIndex[index]+dels.delLen[index] >= self.startPosition + self.seqLen - self.read2Len and dels.delBufferIndex[index] < self.startPosition + self.seqLen)):
                dels.totalReadsCount[index] = dels.totalReadsCount[index] + 1
                if(p < dels.delFraction[index]+0.01):
                    delIndex = dels.delBufferIndex[index]-self.startPosition
                    self.seq[delIndex:delIndex+dels.delLen[index]] = '-' * dels.delLen[index]
                    dels.addedCount[index] = dels.addedCount[index]+1
                    self.suffix = self.suffix + '_' + str(dels.absolutePosition[index])
        self.seq = filter(lambda x: x !='-',self.seq)
        
    def addInsertion(self):
        pass
    
    def pairedEndConverter(self):
        self.read1 = self.seq[0:self.read1Len]
        baseComplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        tmpRead2 = self.seq[-self.read2Len:]
        self.read2 = [baseComplement[x] for x in tmpRead2]
        self.read2.reverse()
        self.read1Score = self.score[0:self.read1Len]
        self.read2Score = self.score[-self.read2Len:]
        self.read2Score.reverse()        

    def writeToFile(self,readsFileHandle):
        readsFileHandle.reads1F.write(self.id + self.suffix + '/1\n')
        readsFileHandle.reads1F.write(''.join(self.read1) + '\n')
        readsFileHandle.reads1F.write('+' + '\n')
        readsFileHandle.reads1F.write(''.join(self.read1Score) + '\n')
        readsFileHandle.reads2F.write(self.id + self.suffix + '/2\n')
        readsFileHandle.reads2F.write(''.join(self.read2) + '\n')
        readsFileHandle.reads2F.write('+' + '\n')
        readsFileHandle.reads2F.write(''.join(self.read2Score) + '\n')
    
    def printer(self):
        print self.id
        print self.seq
        print self.score


class SNP:
    def __init__(self,snpPercentage,sequence,absoluteStart,strand,recordID,snpFraction,bufferRegion):
        #random.seed(10) #delete later
        self.sequence = sequence[bufferRegion:-bufferRegion]
        self.snpCount = int(math.floor(len(self.sequence) * snpPercentage))
        self.snpIndex = sorted(random.sample(range(len(self.sequence)),self.snpCount))
        self.snpBufferIndex = [x + bufferRegion for x in self.snpIndex]
        self.snpAlleles = self.findAllele(self.snpIndex,self.sequence)
        self.snpFraction = [snpFraction] * self.snpCount    
        self.addedCount = [0] * len(self.snpIndex)
        self.recordID = recordID
        self.strand = strand
        self.totalReadsCount = [0] * len(self.snpIndex)
        self.absolutePosition = [x + absoluteStart for x in self.snpBufferIndex]

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

class Deletion:
    def __init__(self,delPercentage,sequence,absoluteStart,strand,recordID,delFraction,bufferRegion,delLen):
        self.sequence = sequence[bufferRegion:-bufferRegion]
        self.delCount = int(math.floor(len(self.sequence) * delPercentage))
        self.delLen = [delLen] * self.delCount
        self.delIndex = sorted(random.sample(range(len(self.sequence)),self.delCount))
        self.delAllele = [self.sequence[self.delIndex[x]:self.delIndex[x]+self.delLen[x]] for x in range(self.delCount)]
        self.delFraction = [delFraction] * self.delCount
        self.delBufferIndex = [x + bufferRegion for x in self.delIndex]
        self.snpFraction = [delFraction] * self.delCount    
        self.addedCount = [0] * len(self.delIndex)
        self.recordID = recordID
        self.strand = strand
        self.totalReadsCount = [0] * len(self.delIndex)
        self.absolutePosition = [x + absoluteStart for x in self.delBufferIndex]
    
    def writeToFile(self,fileH):
        for i in range(len(self.addedCount)):
            fileH.write('\t'.join([self.recordID,str(self.absolutePosition[i]),str(self.delAllele[i]),'.',str(self.addedCount[i]),str(self.totalReadsCount[i])])+'\n')

    def printer(self):
        print self.delIndex
        print self.delAlleles
        print self.delFraction	

    
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
    def __init__(self,filePrefix,snpOutFile,errorOutFile):
        self.fileName = [filePrefix+'1.fastq',filePrefix+'2.fastq']
        self.snpOutFile = snpOutFile
        self.errorOutFile = errorOutFile

    def fileOpen(self):
        self.errorFH = open(self.errorOutFile,'w')
        self.snpFH = open(self.snpOutFile,'w')
        self.reads1F = open(self.fileName[0],'a')
        self.reads2F = open(self.fileName[1],'a')
            
    def fileClose(self):
        self.errorFH.close()
        self.snpFH.close()
        self.reads1F.close()
        self.reads2F.close()
    
    def fileDelete(self):
        for file in self.fileName:
            if os.path.isfile(file):
                print 'Removing ' + str(file) + '!'
                os.system('rm ' + file)
        if os.path.isfile(self.snpOutFile):
            print 'Removing ' + self.snpOutFile + '!'
            os.system('rm ' + self.snpOutFile)

        if os.path.isfile(self.errorOutFile):
            print 'Removing ' + self.errorOutFile + '!'
            os.system('rm ' + self.errorOutFile)


# Function for generating reads
def generateReads(seqBiopython,readsFileH,parameters):
    sequence = seqBiopython.seq
    idPrefix = seqBiopython.id
    geneStart = int(idPrefix.split('_')[3])
    strand = idPrefix.split('_')[1]
    #readLen = parameters['readLen']
    bufferRegion = parameters['bufferRegion']
    snp = SNP(parameters['snpPercentage'],sequence,geneStart,strand,idPrefix,parameters['snpFraction'],bufferRegion)
    dels = Deletion(parameters['snpPercentage'],sequence,geneStart,strand,idPrefix,parameters['snpFraction'],bufferRegion,parameters['delLen'])
    error = Error()
    seqLen =len(sequence)
#    for i in range(seqLen - min(readLen,read1Len+read2Len)):
    for i in range(0,seqLen-bufferRegion):
        read1Len = parameters['read1Len']
        read2Len = parameters['read2Len']
        multiple = parameters['coverage']/float(read1Len+read2Len)
        while(multiple > 0):
            p = random.uniform(0,1)
            if(p <= multiple):
                idPrefixSuffix = idPrefix + '_' + str(multiple)            
	        #dynamicGap = random.sample(range(min(readLen,seqLen - i)-read1Len-read2Len),1)[0]
                read = Read(i,sequence,idPrefixSuffix,error,parameters['scoreRange'],parameters['errorRate'],read1Len,read2Len)
                #read.addSNP(snp)
                read.addDeletion(dels)
                read.pairedEndConverter()
                read.writeToFile(readsFileH)
            multiple = multiple-1
    #snp.writeToFile(readsFileH.snpFH)
    dels.writeToFile(readsFileH.snpFH)
    error.writeToFile(readsFileH.errorFH)

def argumentParser(argv):
    try:
        opts, args = getopt.getopt(argv,'hG:S:F:E:12',['help','coverage=','readsOutFile=','errorRate=','readLength=','snpFraction=','scoreRange=','bufferRegion=','snpPercentage=','snpOutFile=','refGenome=','readsOutFile=','errorOutFile=','delLength='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    parameters = {'refGenome' : '',
        'snpOutFile' : '',
        'errorOutFile' : '',
        'readsOutFilePrefix' : '',
        'coverage' : 50,
        'errorRate' : 0.0,
        'readLen' : 400,
        'snpPercentage' : 0.01,
        'snpFraction' : 0.5,
        'scoreRange' : [0,40],
        'bufferRegion' : 100,
        'read1Len' : 100,
        'read2Len' : 100,
        'delLen' : 1,
    }
	
    for opt,arg in opts:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ('-G','--refGenome'):
            parameters['refGenome'] = arg
        elif opt in ('-S','--snpOutFile'):
            parameters['snpOutFile'] = arg
        elif opt in ('-F','--readsOutFile'):
            parameters['readsOutFilePrefix'] = arg
        elif opt in ('-E','--errorOutFile'):
            parameters['errorOutFile'] = arg
        elif opt == '--errorRate':
            parameters['errorRate'] = float(arg)
        elif opt == '--readLength':
            parameters['readLen'] = int(arg)
        elif opt == '--snpFraction':
            parameters['snpFraction'] = float(arg)
        elif opt == '--scoreRange':
            parameters['scoreRange'] = arg
        elif opt == '--bufferRegion':
            parameters['bufferRegion'] = int(arg)
        elif opt == '--snpPercentage':
            parameters['snpPercentage'] = float(arg)
        elif opt == '--coverage':
            parameters['coverage'] = int(arg)
        elif opt == '-1':
            parameters['read1Len'] = int(arg)
        elif opt == '-2':
            parameters['read2Len'] = int(arg)
        elif opt =='--delLength':
            parameters['delLen'] = int(arg)
        else:
            assert False, 'unhandeld option'
    
    return parameters
			
def usage():
    pass
    

def main():
    para = argumentParser(sys.argv[1:])
    #print para
    readsFileH = ReadsFileHandle(para['readsOutFilePrefix'],para['snpOutFile'],para['errorOutFile'])
    readsFileH.fileDelete()
    readsFileH.fileOpen()
    print 'Start simulating.'
    for seqRecord in SeqIO.parse(para['refGenome'],'fasta'):
        generateReads(seqRecord,readsFileH,para)
    readsFileH.fileClose()

if __name__ == "__main__":
    main()
