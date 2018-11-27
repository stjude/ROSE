#130428

#ROSE_geneMapper.py

#main method wrapped script to take the enhancer region table output of ROSE_Main and map genes to it
#will create two outputs a gene mapped region table where each row is an enhancer
#and a gene table where each row is a gene
#does this by default for super-enhancers only

import sys



import ROSE_utils


import os

from string import upper,join

from collections import defaultdict


#==================================================================
#====================MAPPING GENES TO ENHANCERS====================
#==================================================================



def mapEnhancerToGene(annotFile,enhancerFile,transcribedFile='',uniqueGenes=True,searchWindow =50000,noFormatTable = False):
    
    '''
    maps genes to enhancers. if uniqueGenes, reduces to gene name only. Otherwise, gives for each refseq
    '''
    startDict = ROSE_utils.makeStartDict(annotFile)
    enhancerTable = ROSE_utils.parseTable(enhancerFile,'\t')

    #internal parameter for debugging
    byRefseq = False


    if len(transcribedFile) > 0:
        transcribedTable = ROSE_utils.parseTable(transcribedFile,'\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = startDict.keys()

    print('MAKING TRANSCRIPT COLLECTION')
    transcribedCollection = ROSE_utils.makeTranscriptCollection(annotFile,0,0,500,transcribedGenes)


    print('MAKING TSS COLLECTION')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,0,0))


    #this turns the tssLoci list into a LocusCollection
    #50 is the internal parameter for LocusCollection and doesn't really matter
    tssCollection = ROSE_utils.LocusCollection(tssLoci,50)

    

    geneDict = {'overlapping':defaultdict(list),'proximal':defaultdict(list)}

    #dictionaries to hold ranks and superstatus of gene nearby enhancers
    rankDict = defaultdict(list)
    superDict= defaultdict(list)

    #list of all genes that appear in this analysis
    overallGeneList = []

    if noFormatTable:
        #set up the output tables
        #first by enhancer
        enhancerToGeneTable = [enhancerTable[0]+['OVERLAP_GENES','PROXIMAL_GENES','CLOSEST_GENE']]

        
    else:
        #set up the output tables
        #first by enhancer
        enhancerToGeneTable = [enhancerTable[0][0:9]+['OVERLAP_GENES','PROXIMAL_GENES','CLOSEST_GENE'] + enhancerTable[5][-2:]]

        #next by gene
        geneToEnhancerTable = [['GENE_NAME','REFSEQ_ID','PROXIMAL_ENHANCERS']]

    #next make the gene to enhancer table
    geneToEnhancerTable = [['GENE_NAME','REFSEQ_ID','PROXIMAL_ENHANCERS','ENHANCER_RANKS','IS_SUPER']]

        


    for line in enhancerTable:
        if line[0][0] =='#' or line[0][0] == 'R':
            continue

        enhancerString = '%s:%s-%s' % (line[1],line[2],line[3])
        
        enhancerLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.',line[0])

        #overlapping genes are transcribed genes whose transcript is directly in the stitchedLocus         
        overlappingLoci = transcribedCollection.getOverlap(enhancerLocus,'both')           
        overlappingGenes =[]
        for overlapLocus in overlappingLoci:                
            overlappingGenes.append(overlapLocus.ID())

        #proximalGenes are transcribed genes where the tss is within 50kb of the boundary of the stitched loci
        proximalLoci = tssCollection.getOverlap(ROSE_utils.makeSearchLocus(enhancerLocus,searchWindow,searchWindow),'both')           
        proximalGenes =[]
        for proxLocus in proximalLoci:
            proximalGenes.append(proxLocus.ID())


        distalLoci = tssCollection.getOverlap(ROSE_utils.makeSearchLocus(enhancerLocus,1000000,1000000),'both')           
        distalGenes =[]
        for proxLocus in distalLoci:
            distalGenes.append(proxLocus.ID())

            
            
        overlappingGenes = ROSE_utils.uniquify(overlappingGenes)
        proximalGenes = ROSE_utils.uniquify(proximalGenes)
        distalGenes = ROSE_utils.uniquify(distalGenes)
        allEnhancerGenes = overlappingGenes + proximalGenes + distalGenes
        #these checks make sure each gene list is unique.
        #technically it is possible for a gene to be overlapping, but not proximal since the
        #gene could be longer than the 50kb window, but we'll let that slide here
        for refID in overlappingGenes:
            if proximalGenes.count(refID) == 1:
                proximalGenes.remove(refID)

        for refID in proximalGenes:
            if distalGenes.count(refID) == 1:
                distalGenes.remove(refID)


        #Now find the closest gene
        if len(allEnhancerGenes) == 0:
            closestGene = ''
        else:
            #get enhancerCenter
            enhancerCenter = (int(line[2]) + int(line[3]))/2

            #get absolute distance to enhancer center
            distList = [abs(enhancerCenter - startDict[geneID]['start'][0]) for geneID in allEnhancerGenes]
            #get the ID and convert to name
            closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']

        #NOW WRITE THE ROW FOR THE ENHANCER TABLE
        if noFormatTable:

            newEnhancerLine = list(line)
            newEnhancerLine.append(join(ROSE_utils.uniquify([startDict[x]['name'] for x in overlappingGenes]),','))
            newEnhancerLine.append(join(ROSE_utils.uniquify([startDict[x]['name'] for x in proximalGenes]),','))
            newEnhancerLine.append(closestGene)

        else:
            newEnhancerLine = line[0:9]
            newEnhancerLine.append(join(ROSE_utils.uniquify([startDict[x]['name'] for x in overlappingGenes]),','))
            newEnhancerLine.append(join(ROSE_utils.uniquify([startDict[x]['name'] for x in proximalGenes]),','))
            newEnhancerLine.append(closestGene)
            newEnhancerLine += line[-2:]

        enhancerToGeneTable.append(newEnhancerLine)
        #Now grab all overlapping and proximal genes for the gene ordered table

        overallGeneList +=overlappingGenes
        for refID in overlappingGenes:
            geneDict['overlapping'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))
            
        overallGeneList+=proximalGenes
        for refID in proximalGenes:
            geneDict['proximal'][refID].append(enhancerString)
            rankDict[refID].append(int(line[-2]))
            superDict[refID].append(int(line[-1]))



    #End loop through
    
    #Make table by gene
    overallGeneList = ROSE_utils.uniquify(overallGeneList)  

    #use enhancer rank to order
    rankOrder = ROSE_utils.order([min(rankDict[x]) for x in overallGeneList])
        
    usedNames = []
    for i in rankOrder:
        refID = overallGeneList[i]
        geneName = startDict[refID]['name']
        if usedNames.count(geneName) > 0 and uniqueGenes == True:

            continue
        else:
            usedNames.append(geneName)
        
        proxEnhancers = geneDict['overlapping'][refID]+geneDict['proximal'][refID]
        
        superStatus = max(superDict[refID])
        enhancerRanks = join([str(x) for x in rankDict[refID]],',')
    
        newLine = [geneName,refID,join(proxEnhancers,','),enhancerRanks,superStatus]
        geneToEnhancerTable.append(newLine)

    #resort enhancerToGeneTable
    if noFormatTable:
        return enhancerToGeneTable,geneToEnhancerTable
    else:
        enhancerOrder = ROSE_utils.order([int(line[-2]) for line in enhancerToGeneTable[1:]])
        sortedTable = [enhancerToGeneTable[0]]
        for i in enhancerOrder:
            sortedTable.append(enhancerToGeneTable[(i+1)])

        return sortedTable,geneToEnhancerTable



#==================================================================
#=========================MAIN METHOD==============================
#==================================================================

def main():
    '''
    main run call
    '''
    debug = False


    from optparse import OptionParser
    usage = "usage: %prog [options] -g [GENOME] -i [INPUT_ENHANCER_FILE]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-i","--i", dest="input",nargs = 1, default=None,
                      help = "Enter a ROSE ranked enhancer or super-enhancer file")
    parser.add_option("-g","--genome", dest="genome",nargs = 1, default=None,
                      help = "Enter the genome build (MM10,MM9,MM8,HG18,HG19,HG38)")

    #optional flags
    parser.add_option("-l","--list", dest="geneList",nargs = 1, default=None,
                      help = "Enter a gene list to filter through")
    parser.add_option("-o","--out", dest="out",nargs = 1, default=None,
                      help = "Enter an output folder. Default will be same folder as input file")
    parser.add_option("-w","--window", dest="window",nargs = 1, default=50000,
                      help = "Enter a search distance for genes. Default is 50,000bp")
    parser.add_option("-f","--format", dest="formatTable",action= "store_true", default=False,
                      help = "If flagged, maintains original formatting of input table")

    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()


    if not options.input or not options.genome:

        parser.print_help()
        exit()

    #GETTING THE INPUT
    enhancerFile = options.input
    window = int(options.window)

    #making the out folder if it doesn't exist
    if options.out:
        outFolder = ROSE_utils.formatFolder(options.out,True)
    else:
        outFolder = join(enhancerFile.split('/')[0:-1],'/') + '/'


    #GETTING THE GENOME
    genome = options.genome
    print('USING %s AS THE GENOME' % genome)

    #CHECK FORMATTING FLAG
    if options.formatTable:
        noFormatTable =True
    else:
        noFormatTable = False

    #GETTING THE CORRECT ANNOT FILE
    cwd = os.getcwd()
    genomeDict = {
        'HG18':'%s/annotation/hg18_refseq.ucsc' % (cwd),
        'MM9': '%s/annotation/mm9_refseq.ucsc' % (cwd),
        'HG19':'%s/annotation/hg19_refseq.ucsc' % (cwd),
	'HG38':'%s/annotation/hg38_refseq.ucsc' % (cwd),
        'MM8': '%s/annotation/mm8_refseq.ucsc' % (cwd),
        'MM10':'%s/annotation/mm10_refseq.ucsc' % (cwd),
        }

    annotFile = genomeDict[upper(genome)]

    #GETTING THE TRANSCRIBED LIST
    if options.geneList:

        transcribedFile = options.geneList
    else:
        transcribedFile = ''

    enhancerToGeneTable,geneToEnhancerTable = mapEnhancerToGene(annotFile,enhancerFile,transcribedFile,True,window,noFormatTable)

    #Writing enhancer output
    enhancerFileName = enhancerFile.split('/')[-1].split('.')[0]

    if window != 50000:
        #writing the enhancer table
        out1 = '%s%s_ENHANCER_TO_GENE_%sKB.txt' % (outFolder,enhancerFileName,window/1000)
        ROSE_utils.unParseTable(enhancerToGeneTable,out1,'\t')

        #writing the gene table
        out2 = '%s%s_GENE_TO_ENHANCER_%sKB.txt' % (outFolder,enhancerFileName,window/1000)
        ROSE_utils.unParseTable(geneToEnhancerTable,out2,'\t')
    else:
        #writing the enhancer table
        out1 = '%s%s_ENHANCER_TO_GENE.txt' % (outFolder,enhancerFileName)
        ROSE_utils.unParseTable(enhancerToGeneTable,out1,'\t')

        #writing the gene table
        out2 = '%s%s_GENE_TO_ENHANCER.txt' % (outFolder,enhancerFileName)
        ROSE_utils.unParseTable(geneToEnhancerTable,out2,'\t')

if __name__ == "__main__":
    main()
