#bamToGFF_turbo.py

#20130716

#script to grab reads from a bam that align to a .gff file

#uses the bamliquidator super fast uber thingy written by Xin Zhou
import sys



from ROSE_utils import *

from collections import defaultdict

import os
import string
import subprocess

        
        



def mapBamToGFF(bamFile,gff,sense = '.',extension = 200,rpm = False,clusterGram = None,matrix = None):
    '''maps reads from a bam to a gff'''

    #creating a new gff to output
    newGFF = []
    #reading in the bam
    bam = Bam(bamFile)

    #getting RPM normalization
    if rpm:    
        MMR= round(float(bam.getTotalReads('mapped'))/1000000,4)
    else:
        MMR = 1

    print('using a MMR value of %s' % (MMR))

    #creating a sense trans 
    senseTrans = maketrans('-+.','+-+')
    
    #reading in the gff
    if type(gff) == str:
        gff = parseTable(gff,'\t')

    #setting up a clustergram table
    if clusterGram:
        binSize = int(clusterGram)
        binSizeList = []
        #now go through each line of the gff and make sure they're all the same length
        for i in range(0,len(gff),1):
            line = gff[i]
            gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
            binSizeList.append(gffLocus.len()/binSize)
        binSizeList = uniquify(binSizeList)
        if len(binSizeList) > 1: 
            print('WARNING: lines in gff are of different length. Output clustergram will have variable row length')
        newGFF.append(['GENE_ID','locusLine'] + [str(x*binSize)+'_'+bamFile.split('/')[-1] for x in range(1,max(binSizeList)+1,1)])        
        
    #setting up a maxtrix table
    if matrix:
        newGFF.append(['GENE_ID','locusLine'] + ['bin_'+str(n)+'_'+bamFile.split('/')[-1] for n in range(1,int(matrix)+1,1)])
        nBin = int(matrix)

    #getting and processing reads for gff lines
    ticker = 0
    print('Number lines processed')
    for line in gff:
        line = line[0:9]
        if ticker%100 == 0:
            print ticker
        ticker+=1
        gffLocus = Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
        
        #get the nBin and binSize
        if clusterGram:
            nBin =gffLocus.len()/int(clusterGram)
            binSize = int(clusterGram)
        if matrix:
            nBin = int(matrix)
            binSize = gffLocus.len()/nBin
            
        #flippy flip if sense is negative
        if sense == '-':
            bamSense = string.translate(gffLocus.sense(),senseTrans)
        elif sense == '+':
            bamSense = gffLocus.sense()
        else:
            bamSense = '.'
        #using the bamLiquidator to get the readstring            
        #print('using nBin of %s' % nBin)
        bamliquidatorString = '/usr/bin/bamliquidator'
        bamCommand = "%s %s %s %s %s %s %s %s" % (bamliquidatorString,bamFile,line[0],gffLocus.start(),gffLocus.end(),bamSense,nBin,extension)
        #print(bamCommand)
        getReads = subprocess.Popen(bamCommand,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)
        readString = getReads.communicate()
        denList = readString[0].split('\n')[:-1]
        
        #flip the denList if the actual gff region is -
        if gffLocus.sense() == '-':
            denList = denList[::-1]

        #converting from units of total bp of read sequence per bin to rpm/bp
        denList = [round(float(x)/binSize/MMR,4) for x in denList]
        
        #if the gff region is - strand, flip the

        clusterLine = [gffLocus.ID(),gffLocus.__str__()] + denList
        newGFF.append(clusterLine)
            
    return newGFF
        
                
                
            


            
                
    
def convertEnrichedRegionsToGFF(enrichedRegionFile):
    '''converts a young lab enriched regions file into a gff'''
    newGFF = []
    enrichedRegions = open(enrichedRegionFile,'r')
    header = enrichedRegions.readline()
    i = 0
    for line in enrichedRegions:
        line = line[:-1].split('\t')
        newLine = ['chr'+line[0],'row_'+str(i),line[4],line[1],line[2],'','.','','row_'+str(i),'']
        newGFF.append(newLine)
        i+=1
    return newGFF

        
#python bamToGFF.py --density --floor 0 -b test.sam.sorted.bam -g pol2_sample.gff -o pol2_sample_mapped.gff

def main():
    from optparse import OptionParser
    usage = "usage: %prog [options] -b [SORTED BAMFILE] -i [INPUTFILE] -o [OUTPUTFILE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-b","--bam", dest="bam",nargs = 1, default=None,
                      help = "Enter .bam file to be processed.")
    parser.add_option("-i","--input", dest="input",nargs = 1, default=None,
                      help = "Enter .gff or ENRICHED REGION file to be processed.")
    #output flag
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter the output filename.")
    #additional options
    parser.add_option("-s","--sense", dest="sense",nargs = 1, default='.',
                      help = "Map to '+','-' or 'both' strands. Default maps to both.")
    parser.add_option("-e","--extension", dest="extension",nargs = 1, default=200,
                      help = "Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=False,
                      help = "Normalizes density to reads per million (rpm)")
    parser.add_option("-c","--cluster", dest="cluster",nargs = 1, default=None,
                      help = "Outputs a fixed bin size clustergram. user must specify bin size.")
    parser.add_option("-m","--matrix", dest="matrix",nargs = 1, default=None,
                      help = "Outputs a variable bin sized matrix. User must specify number of bins.")
    (options,args) = parser.parse_args()

    print(options)
    print(args)

   
    if options.sense:
        if ['+','-','.','both'].count(options.sense) == 0:
            print('ERROR: sense flag must be followed by +,-,.,both')
            parser.print_help()
            exit()

    if options.cluster and options.matrix:
        print('ERROR: Cannot specify both matrix and clustergram flags.')
        parser.print_help()
        exit()

    if options.matrix:
        try:
            int(options.matrix)
        except:
            print('ERROR: User must specify an integer bin number for matrix (try 50)')
            parser.print_help()
            exit()
            
    if options.cluster:
        try:
            int(options.cluster)
        except:
            print('ERROR: User must specify an integer bin size for clustergram (try 25)')
            parser.print_help()
            exit()

    
    
    if options.input and options.bam:
        inputFile = options.input
        if inputFile.split('.')[-1] != 'gff':
            print('converting file to a .gff')
            gffFile = convertEnrichedRegionsToGFF(inputFile)
        else:
            gffFile = inputFile

        bamFile = options.bam
        
        if options.output == None:
            output = os.getcwd() + inputFile.split('/')[-1]+'.mapped'
        else:
            output = options.output
        if options.cluster:
            print('mapping to GFF and making clustergram with fixed bin width')
            newGFF = mapBamToGFF(bamFile,gffFile,options.sense,int(options.extension),options.rpm,int(options.cluster),None)
        elif options.matrix:
            print('mapping to GFF and making a matrix with fixed bin number')
            newGFF = mapBamToGFF(bamFile,gffFile,options.sense,int(options.extension),options.rpm,None,int(options.matrix))
            
        unParseTable(newGFF,output,'\t')
    else:
        parser.print_help()
        

        
if __name__ == "__main__":
    main()
