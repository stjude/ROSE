
#bamToGFF.py

#script to grab reads from a bam that align to a .gff file
import sys
import re

import ROSE_utils


from collections import defaultdict

import os

from string import join,upper,maketrans



#=====================================================================
#====================MAPPING BAM READS TO GFF REGIONS=================
#=====================================================================


def mapBamToGFF(bamFile,gff,sense = 'both',extension = 200,floor = 0,rpm = False,matrix = None):

#def mapBamToGFF(bamFile,gff,sense = 'both',unique = 0,extension = 200,floor = 0,density = False,rpm = False,binSize = 25,clusterGram = None,matrix = None,raw = False,includeJxnReads = False):
    '''maps reads from a bam to a gff'''
    floor = int(floor)
    
    #USING BAM CLASS
    bam = ROSE_utils.Bam(bamFile)


    #new GFF to write to
    newGFF = []
    #millionMappedReads


    if rpm:    
        MMR= round(float(bam.getTotalReads('mapped'))/1000000,4)
    else:
        MMR = 1

    print('using a MMR value of %s' % (MMR))
    
    senseTrans = maketrans('-+.','+-+')

    if ROSE_utils.checkChrStatus(bamFile) == 1:
      print "has chr"
      hasChrFlag = 1
      #sys.exit();
    else:
      print "does not have chr"
      hasChrFlag = 0
      #sys.exit()
      
    if type(gff) == str:
        gff = ROSE_utils.parseTable(gff,'\t')
        
    #setting up a maxtrix table

    newGFF.append(['GENE_ID','locusLine'] + ['bin_'+str(n)+'_'+bamFile.split('/')[-1] for n in range(1,int(matrix)+1,1)])        

    #getting and processing reads for gff lines
    ticker = 0
    print('Number lines processed')
    for line in gff:
        line = line[0:9]
        if ticker%100 == 0:
            print ticker
        ticker+=1
        if not hasChrFlag:
	  line[0] = re.sub(r"chr",r"",line[0])
        gffLocus = ROSE_utils.Locus(line[0],int(line[3]),int(line[4]),line[6],line[1])
        #print line[0]
        #sys.exit()
        searchLocus = ROSE_utils.makeSearchLocus(gffLocus,int(extension),int(extension))
        
        reads = bam.getReadsLocus(searchLocus,'both',False,'none')
        #now extend the reads and make a list of extended reads
        extendedReads = []
        for locus in reads:
            if locus.sense() == '+' or locus.sense() == '.':
                locus = ROSE_utils.Locus(locus.chr(),locus.start(),locus.end()+extension,locus.sense(), locus.ID())
            if locus.sense() == '-':
                locus = ROSE_utils.Locus(locus.chr(),locus.start()-extension,locus.end(),locus.sense(),locus.ID())
            extendedReads.append(locus)
        if gffLocus.sense() == '+' or gffLocus.sense == '.':
            senseReads = filter(lambda x:x.sense() == '+' or x.sense() == '.',extendedReads)
            antiReads = filter(lambda x:x.sense() == '-',extendedReads)
        else:
            senseReads = filter(lambda x:x.sense() == '-' or x.sense() == '.',extendedReads)
            antiReads = filter(lambda x:x.sense() == '+',extendedReads)

        senseHash = defaultdict(int)
        antiHash = defaultdict(int)

        #filling in the readHashes             
        if sense == '+' or sense == 'both' or sense =='.':
            for read in senseReads:
                for x in range(read.start(),read.end()+1,1):
                    senseHash[x]+=1
        if sense == '-' or sense == 'both' or sense == '.':
            #print('foo')
            for read in antiReads:
                for x in range(read.start(),read.end()+1,1):
                    antiHash[x]+=1

        #now apply flooring and filtering for coordinates
        keys = ROSE_utils.uniquify(senseHash.keys()+antiHash.keys())
        if floor > 0:

            keys = filter(lambda x: (senseHash[x]+antiHash[x]) > floor,keys)
        #coordinate filtering
        keys = filter(lambda x: gffLocus.start() < x < gffLocus.end(),keys)


        #setting up the output table
        clusterLine = [gffLocus.ID(),gffLocus.__str__()]

        #getting the binsize
        binSize = (gffLocus.len()-1)/int(matrix)
        nBins = int(matrix)
        if binSize == 0:
            clusterLine+=['NA']*int(matrix)
            newGFF.append(clusterLine)
            continue
        n=0
        if gffLocus.sense() == '+' or gffLocus.sense() =='.' or gffLocus.sense() == 'both':
            i = gffLocus.start()

            while n <nBins:
                n+=1
                binKeys = filter(lambda x: i < x < i+binSize,keys)
                binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys]))/binSize
                clusterLine+=[round(binDen/MMR,4)]
                i = i+binSize
        else:
            i = gffLocus.end()
            while n < nBins:
                n+=1
                binKeys = filter(lambda x: i-binSize < x < i,keys)
                binDen = float(sum([senseHash[x]+antiHash[x] for x in binKeys]))/binSize
                clusterLine+=[round(binDen/MMR,4)]
                i = i-binSize
        newGFF.append(clusterLine)
        
            
    return newGFF
        
                
                
#=====================================================================
#============================MAIN METHOD==============================
#=====================================================================
        

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
    parser.add_option("-s","--sense", dest="sense",nargs = 1, default='both',
                      help = "Map to '+','-' or 'both' strands. Default maps to both.")


    parser.add_option("-f","--floor", dest="floor",nargs =1, default=0,
                      help = "Sets a read floor threshold necessary to count towards density")    
    parser.add_option("-e","--extension", dest="extension",nargs = 1, default=200,
                      help = "Extends reads by n bp. Default value is 200bp")
    parser.add_option("-r","--rpm", dest="rpm",action = 'store_true', default=False,
                      help = "Normalizes density to reads per million (rpm)")


    parser.add_option("-m","--matrix", dest="matrix",nargs = 1, default=None,
                      help = "Outputs a variable bin sized matrix. User must specify number of bins.")

    (options,args) = parser.parse_args()

    print(options)
    print(args)

    if options.bam:
        bamFile = options.bam
        fullPath = os.path.abspath(bamFile)
        bamName = fullPath.split('/')[-1].split('.')[0]
        pathFolder = join(fullPath.split('/')[0:-1],'/')
        fileList = os.listdir(pathFolder)
        hasBai = False
        for fileName in fileList:
            if fileName.count(bamName) == 1 and fileName.count('.bai') == 1:
                hasBai = True

        if not hasBai:
            print('ERROR: no associated .bai file found with bam. Must use a sorted bam with accompanying index file')
            parser.print_help()
            exit()
   
    if options.sense:
        if ['+','-','.','both'].count(options.sense) == 0:
            print('ERROR: sense flag must be followed by +,-,.,both')
            parser.print_help()
            exit()


    if options.matrix:
        try:
            int(options.matrix)
        except:
            print('ERROR: User must specify an integer bin number for matrix (try 50)')
            parser.print_help()
            exit()
            

    
    
    if options.input and options.bam:
        inputFile = options.input
        gffFile = inputFile

        bamFile = options.bam
        
        if options.output == None:
            output = os.getcwd() + inputFile.split('/')[-1]+'.mapped'
        else:
            output = options.output
        if options.matrix:
            print('mapping to GFF and making a matrix with fixed bin number')

            newGFF = mapBamToGFF(bamFile,gffFile,options.sense,int(options.extension),options.floor,options.rpm,options.matrix)

            
        ROSE_utils.unParseTable(newGFF,output,'\t')
    else:
        parser.print_help()
                
 
if __name__ == "__main__":
    main()
