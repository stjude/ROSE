#mapEnhancerFromFactor.py
'''
PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS,
AND RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS
APRIL 11, 2013
VERSION 0.1
CONTACT: youngcomputation@wi.mit.edu
'''

import sys



import ROSE_utils

import time

import os

from string import upper,join

from collections import defaultdict

#==================================================================
#=====================REGION STITCHING=============================
#==================================================================


def regionStitching(inputGFF,stitchWindow,tssWindow,annotFile,removeTSS=True):
    print('PERFORMING REGION STITCHING')
    #first have to turn bound region file into a locus collection

    #need to make sure this names correctly... each region should have a unique name
    boundCollection = ROSE_utils.gffToLocusCollection(inputGFF)

    debugOutput = []
    #filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if removeTSS:
        #first make a locus collection of TSS
        startDict = ROSE_utils.makeStartDict(annotFile)

        #now makeTSS loci for active genes
        removeTicker=0
        #this loop makes a locus centered around +/- tssWindow of transcribed genes
        #then adds it to the list tssLoci
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,tssWindow,tssWindow))


        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)

        #gives all the loci in boundCollection
        boundLoci = boundCollection.getLoci()

        #this loop will check if each bound region is contained by the TSS exclusion zone
        #this will drop out a lot of the promoter only regions that are tiny
        #typical exclusion window is around 2kb
        for locus in boundLoci:
            if len(tssCollection.getContainers(locus,'both'))>0:
                
                #if true, the bound locus overlaps an active gene
                boundCollection.remove(locus)
                debugOutput.append([locus.__str__(),locus.ID(),'CONTAINED'])
                removeTicker+=1
        print('REMOVED %s LOCI BECAUSE THEY WERE CONTAINED BY A TSS' % (removeTicker))

    #boundCollection is now all enriched region loci that don't overlap an active TSS
    stitchedCollection = boundCollection.stitchCollection(stitchWindow,'both')

    if removeTSS:
        #now replace any stitched region that overlap 2 distinct genes
        #with the original loci that were there
        fixedLoci = []
        tssLoci = []
        for geneID in startDict.keys():
            tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,50,50))


        #this turns the tssLoci list into a LocusCollection
        #50 is the internal parameter for LocusCollection and doesn't really matter
        tssCollection = ROSE_utils.LocusCollection(tssLoci,50)
        removeTicker = 0
        originalTicker = 0
        for stitchedLocus in stitchedCollection.getLoci():
            overlappingTSSLoci = tssCollection.getOverlap(stitchedLocus,'both')
            tssNames = [startDict[tssLocus.ID()]['name'] for tssLocus in overlappingTSSLoci]
            tssNames = ROSE_utils.uniquify(tssNames)
            if len(tssNames) > 2:
            
                #stitchedCollection.remove(stitchedLocus)
                originalLoci = boundCollection.getOverlap(stitchedLocus,'both')
                originalTicker+=len(originalLoci)
                fixedLoci+=originalLoci
                debugOutput.append([stitchedLocus.__str__(),stitchedLocus.ID(),'MULTIPLE_TSS'])
                removeTicker+=1
            else:
                fixedLoci.append(stitchedLocus)

        print('REMOVED %s STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs' % (removeTicker))
        print('ADDED BACK %s ORIGINAL LOCI' % (originalTicker))
        fixedCollection = ROSE_utils.LocusCollection(fixedLoci,50)
        return fixedCollection,debugOutput
    else:
        return stitchedCollection,debugOutput

#==================================================================
#=====================REGION LINKING MAPPING=======================
#==================================================================

def mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,output,refName):


    '''
    makes a table of factor density in a stitched locus and ranks table by number of loci stitched together
    '''

    
    print('FORMATTING TABLE')
    loci = stitchedCollection.getLoci()

    locusTable = [['REGION_ID','CHROM','START','STOP','NUM_LOCI','CONSTITUENT_SIZE']]

    lociLenList = []

    #strip out any that are in chrY
    for locus in list(loci):
        if locus.chr() == 'chrY':
            loci.remove(locus)
    
    for locus in loci:
        #numLociList.append(int(stitchLocus.ID().split('_')[1]))
        lociLenList.append(locus.len())
        #numOrder = order(numLociList,decreasing=True)
    lenOrder = ROSE_utils.order(lociLenList,decreasing=True)
    ticker = 0
    for i in lenOrder:
        ticker+=1
        if ticker%1000 ==0:
            print(ticker)
        locus = loci[i]

        #First get the size of the enriched regions within the stitched locus
        refEnrichSize = 0
        refOverlappingLoci = referenceCollection.getOverlap(locus,'both')
        for refLocus in refOverlappingLoci:
            refEnrichSize+=refLocus.len()

        try:
            stitchCount = int(locus.ID().split('_')[0])
        except ValueError:
            stitchCount = 1
        
        locusTable.append([locus.ID(),locus.chr(),locus.start(),locus.end(),stitchCount,refEnrichSize])
        
            

    print('GETTING MAPPED DATA')
    for bamFile in bamFileList:
        
        bamFileName = bamFile.split('/')[-1]

        print('GETTING MAPPING DATA FOR  %s' % bamFile)
        #assumes standard convention for naming enriched region gffs
        
        #opening up the mapped GFF
        print('OPENING %s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName))

        mappedGFF =ROSE_utils.parseTable('%s%s_%s_MAPPED.gff' % (mappedFolder,refName,bamFileName),'\t')        

        signalDict = defaultdict(float)
        print('MAKING SIGNAL DICT FOR %s' % (bamFile))
        mappedLoci = []
        for line in mappedGFF[1:]:

            chrom = line[1].split('(')[0]
            start = int(line[1].split(':')[-1].split('-')[0])
            end = int(line[1].split(':')[-1].split('-')[1])
            mappedLoci.append(ROSE_utils.Locus(chrom,start,end,'.',line[0]))
            try:
                signalDict[line[0]] = float(line[2])*(abs(end-start))
            except ValueError:
                print('WARNING NO SIGNAL FOR LINE:')
                print(line)
                continue
                
                
        
        mappedCollection = ROSE_utils.LocusCollection(mappedLoci,500)
        locusTable[0].append(bamFileName)

        for i in range(1,len(locusTable)):
            signal=0.0
            line = locusTable[i]
            lineLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.')
            overlappingRegions = mappedCollection.getOverlap(lineLocus,sense='both')
            for region in overlappingRegions:
                signal+= signalDict[region.ID()]
            locusTable[i].append(signal)

    ROSE_utils.unParseTable(locusTable,output,'\t')



#==================================================================
#=========================MAIN METHOD==============================
#==================================================================

def main():
    '''
    main run call
    '''
    debug = False


    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-i","--i", dest="input",nargs = 1, default=None,
                      help = "Enter a .gff or .bed file of binding sites used to make enhancers")
    parser.add_option("-r","--rankby", dest="rankby",nargs = 1, default=None,
                      help = "bamfile to rank enhancer by")
    parser.add_option("-o","--out", dest="out",nargs = 1, default=None,
                      help = "Enter an output folder")
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (MM9,MM8,HG18,HG19)")
    
    #optional flags
    parser.add_option("-b","--bams", dest="bams",nargs = 1, default=None,
                      help = "Enter a comma separated list of additional bam files to map to")
    parser.add_option("-c","--control", dest="control",nargs = 1, default=None,
                      help = "bamfile to rank enhancer by")
    parser.add_option("-s","--stitch", dest="stitch",nargs = 1, default=12500,
                      help = "Enter a max linking distance for stitching")
    parser.add_option("-t","--tss", dest="tss",nargs = 1, default=0,
                      help = "Enter a distance from TSS to exclude. 0 = no TSS exclusion")




    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()


    if not options.input or not options.rankby or not options.out or not options.genome:
        print('hi there')
        parser.print_help()
        exit()

    #making the out folder if it doesn't exist
    outFolder = ROSE_utils.formatFolder(options.out,True)

    
    #figuring out folder schema
    gffFolder = ROSE_utils.formatFolder(outFolder+'gff/',True)
    mappedFolder = ROSE_utils.formatFolder(outFolder+ 'mappedGFF/',True)


    #GETTING INPUT FILE
    if options.input.split('.')[-1] == 'bed':
        #CONVERTING A BED TO GFF
        inputGFFName = options.input.split('/')[-1][0:-4]
        inputGFFFile = '%s%s.gff' % (gffFolder,inputGFFName)
        ROSE_utils.bedToGFF(options.input,inputGFFFile)
    elif options.input.split('.')[-1] =='gff':
        #COPY THE INPUT GFF TO THE GFF FOLDER
	inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile,gffFolder))        

    else:
        print('WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT')
        #COPY THE INPUT GFF TO THE GFF FOLDER
	inputGFFFile = options.input
        os.system('cp %s %s' % (inputGFFFile,gffFolder))        


    #GETTING THE LIST OF BAMFILES TO PROCESS
    if options.control:        
        bamFileList = [options.rankby,options.control]

    else:
        bamFileList = [options.rankby]

    if options.bams:
        bamFileList += options.bams.split(',')
        bamFileLIst = ROSE_utils.uniquify(bamFileList)
    #optional args

    #Stitch parameter
    stitchWindow = int(options.stitch)
    
    #tss options
    tssWindow = int(options.tss)
    if tssWindow != 0:
        removeTSS = True
    else:
        removeTSS = False

    #GETTING THE BOUND REGION FILE USED TO DEFINE ENHANCERS
    print('USING %s AS THE INPUT GFF' % (inputGFFFile))
    inputName = inputGFFFile.split('/')[-1].split('.')[0]


    #GETTING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)


    #GETTING THE CORRECT ANNOT FILE
    cwd = os.getcwd()
    genomeDict = {
        'HG18':'%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19':'%s/annotation/hg19_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10':'%s/annotation/mm10_refseq.ucsc' % (cwd),
        }

    annotFile = genomeDict[upper(genome)]

    #MAKING THE START DICT
    print('MAKING START DICT')
    startDict = ROSE_utils.makeStartDict(annotFile)


    #LOADING IN THE BOUND REGION REFERENCE COLLECTION
    print('LOADING IN GFF REGIONS')
    referenceCollection = ROSE_utils.gffToLocusCollection(inputGFFFile)


    #NOW STITCH REGIONS
    print('STITCHING REGIONS TOGETHER')
    stitchedCollection,debugOutput = regionStitching(inputGFFFile,stitchWindow,tssWindow,annotFile,removeTSS)

    
    #NOW MAKE A STITCHED COLLECTION GFF
    print('MAKING GFF FROM STITCHED COLLECTION')
    stitchedGFF=ROSE_utils.locusCollectionToGFF(stitchedCollection)
    
    if not removeTSS:
        stitchedGFFFile = '%s%s_%sKB_STITCHED.gff' % (gffFolder,inputName,stitchWindow/1000)
        stitchedGFFName = '%s_%sKB_STITCHED' % (inputName,stitchWindow/1000)
        debugOutFile = '%s%s_%sKB_STITCHED.debug' % (gffFolder,inputName,stitchWindow/1000)
    else:
        stitchedGFFFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.gff' % (gffFolder,inputName,stitchWindow/1000)
        stitchedGFFName = '%s_%sKB_STITCHED_TSS_DISTAL' % (inputName,stitchWindow/1000)
        debugOutFile = '%s%s_%sKB_STITCHED_TSS_DISTAL.debug' % (gffFolder,inputName,stitchWindow/1000)

    #WRITING DEBUG OUTPUT TO DISK
        
    if debug:
        print('WRITING DEBUG OUTPUT TO DISK AS %s' % (debugOutFile))
        ROSE_utils.unParseTable(debugOutput,debugOutFile,'\t')

    #WRITE THE GFF TO DISK
    print('WRITING STITCHED GFF TO DISK AS %s' % (stitchedGFFFile))
    ROSE_utils.unParseTable(stitchedGFF,stitchedGFFFile,'\t')



    #SETTING UP THE OVERALL OUTPUT FILE
    outputFile1 = outFolder + stitchedGFFName + '_ENHANCER_REGION_MAP.txt'

    print('OUTPUT WILL BE WRITTEN TO  %s' % (outputFile1))
    
    #MAPPING TO THE NON STITCHED (ORIGINAL GFF)
    #MAPPING TO THE STITCHED GFF


    # bin for bam mapping
    nBin =1

    #IMPORTANT
    #CHANGE cmd1 and cmd2 TO PARALLELIZE OUTPUT FOR BATCH SUBMISSION
    #e.g. if using LSF cmd1 = "bsub python bamToGFF.py -f 1 -e 200 -r -m %s -b %s -i %s -o %s" % (nBin,bamFile,stitchedGFFFile,mappedOut1)

    for bamFile in bamFileList:

        bamFileName = bamFile.split('/')[-1]

        #MAPPING TO THE STITCHED GFF
        mappedOut1 ='%s%s_%s_MAPPED.gff' % (mappedFolder,stitchedGFFName,bamFileName)
        #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
        cmd1 = "python ROSE_bamToGFF_turbo.py -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,bamFile,stitchedGFFFile,mappedOut1)
        print(cmd1)
        os.system(cmd1)

        #MAPPING TO THE ORIGINAL GFF
        mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
        #WILL TRY TO RUN AS A BACKGROUND PROCESS. BATCH SUBMIT THIS LINE TO IMPROVE SPEED
        cmd2 = "python ROSE_bamToGFF_turbo.py 1 -e 200 -r -m %s -b %s -i %s -o %s &" % (nBin,bamFile,inputGFFFile,mappedOut2)
        print(cmd2)
        os.system(cmd2)
        

    
    print('PAUSING TO MAP')
    time.sleep(10)

    #CHECK FOR MAPPING OUTPUT
    outputDone = False
    ticker = 0
    print('WAITING FOR MAPPING TO COMPLETE. ELAPSED TIME (MIN):')
    while not outputDone:

        '''
        check every 1 minutes for completed output
        '''
        outputDone = True
        if ticker%6 == 0:
            print(ticker*5)
        ticker +=1
        #CHANGE THIS PARAMETER TO ALLOW MORE TIME TO MAP
        if ticker == 120:
            print('ERROR: OPERATION TIME OUT. MAPPING OUTPUT NOT DETECTED AFTER 2 HOURS')
            exit()
            break
        for bamFile in bamFileList:
            
            #GET THE MAPPED OUTPUT NAMES HERE FROM MAPPING OF EACH BAMFILE
            bamFileName = bamFile.split('/')[-1]
            mappedOut1 ='%s%s_%s_MAPPED.gff' % (mappedFolder,stitchedGFFName,bamFileName)

            try:
                 mapFile = open(mappedOut1,'r')
                 mapFile.close()
            except IOError:
                outputDone = False

            mappedOut2 ='%s%s_%s_MAPPED.gff' % (mappedFolder,inputName,bamFileName)
            
            try:
                mapFile = open(mappedOut2,'r')
                mapFile.close()
            except IOError:
                outputDone = False
        if outputDone == True:
            break
        time.sleep(60)
    print('MAPPING TOOK %s MINUTES' % (ticker))

    print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    #CALCULATE DENSITY BY REGION
    mapCollection(stitchedCollection,referenceCollection,bamFileList,mappedFolder,outputFile1,refName = stitchedGFFName)


    time.sleep(10)

    print('CALLING AND PLOTTING SUPER-ENHANCERS')


    if options.control:

        rankbyName = options.rankby.split('/')[-1]
        controlName = options.control.split('/')[-1]
        cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)

    else:
        rankbyName = options.rankby.split('/')[-1]
        controlName = 'NONE'
        cmd = 'R --no-save %s %s %s %s < ROSE_callSuper.R' % (outFolder,outputFile1,inputName,controlName)
    print(cmd)
    os.system(cmd)


    #calling the gene mapper                                                                        
    time.sleep(20)
    superTableFile = "%s_SuperEnhancers.table.txt" % (inputName)
    cmd = "python ROSE_geneMapper.py -g %s -i %s%s" % (genome,outFolder,superTableFile)
    os.system(cmd)    




if __name__ == "__main__":
    main()
