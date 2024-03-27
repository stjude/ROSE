#!/usr/bin/env python3
#130428

#ROSE_geneMapper.py

#main method wrapped script to take the enhancer region table output of ROSE_Main and map genes to it
#will create two outputs a gene mapped region table where each row is an enhancer
#and a gene table where each row is a gene
#does this by default for super-enhancers only

import sys

import ROSE_utils

import os

import string

from collections import defaultdict


#==================================================================
#====================MAPPING GENES TO ENHANCERS====================
#==================================================================



def mapEnhancerToGene(annotFile,enhancerFile,transcribedFile='',uniqueGenes=True,byRefseq=False,subtractInput=False):
    
    '''
    maps genes to enhancers. if uniqueGenes, reduces to gene name only. Otherwise, gives for each refseq
    '''
    print("Herp")
    startDict = ROSE_utils.makeStartDict(annotFile)
    print("Derp")
    enhancerTable = ROSE_utils.parseTable(enhancerFile,'\t')





    if len(transcribedFile) > 0:
        transcribedTable = ROSE_utils.parseTable(transcribedFile,'\t')
        transcribedGenes = [line[1] for line in transcribedTable]
    else:
        transcribedGenes = list(startDict.keys())

    print('MAKING TRANSCRIPT COLLECTION')
    transcribedCollection = ROSE_utils.makeTranscriptCollection(annotFile,0,0,500,transcribedGenes)


    print('MAKING TSS COLLECTION')
    tssLoci = []
    for geneID in transcribedGenes:
        tssLoci.append(ROSE_utils.makeTSSLocus(geneID,startDict,0,0))


    #this turns the tssLoci list into a LocusCollection
    #50 is the internal parameter for LocusCollection and doesn't really matter
    tssCollection = ROSE_utils.LocusCollection(tssLoci,50)



    geneDict = {'overlapping':defaultdict(list),'proximal':defaultdict(list),'enhancerString':defaultdict(list)}
    #list of all genes that appear in this analysis
    overallGeneList = []

    #set up the output tables
    #first by enhancer
    enhancerToGeneTable = [enhancerTable[5][0:6]+['OVERLAP_GENES','PROXIMAL_GENES','CLOSEST_GENE'] + enhancerTable[5][-2:]]
    
    #next by gene
    geneToEnhancerTable = [['GENE_NAME','REFSEQ_ID','PROXIMAL_STITCHED_PEAKS']]

    #have all information
    signalWithGenes = [['GENE_NAME', 'REFSEQ_ID','PROXIMAL_STITCHED_PEAKS', 'SIGNAL']]

    for line in enhancerTable[6:]:

        enhancerString = '%s:%s-%s' % (line[1],line[2],line[3])
        enhancerSignal = int(float(line[6]))
        if subtractInput: enhancerSignal = int(float(line[6]) - float(line[7]))

        enhancerLocus = ROSE_utils.Locus(line[1],line[2],line[3],'.',line[0])


        #overlapping genes are transcribed genes whose transcript is directly in the stitchedLocus         
        overlappingLoci = transcribedCollection.getOverlap(enhancerLocus,'both')           
        overlappingGenes =[]
        for overlapLocus in overlappingLoci:                
            overlappingGenes.append(overlapLocus.ID())

        #proximalGenes are transcribed genes where the tss is within 50kb of the boundary of the stitched loci
        proximalLoci = tssCollection.getOverlap(ROSE_utils.makeSearchLocus(enhancerLocus,50000,50000),'both')           
        proximalGenes =[]
        for proxLocus in proximalLoci:
            proximalGenes.append(proxLocus.ID())


        distalLoci = tssCollection.getOverlap(ROSE_utils.makeSearchLocus(enhancerLocus,50000000,50000000),'both')           
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
            #print enhancerCenter - startDict[geneID]['start'][0]
            #print distList.index(min(distList))
            #print min(distList)
            #print len(distList)
	        #print len(allEnhancerGenes[distList.index(min(distList))])
	        #print line
	        #print len(startDict[allEnhancerGenes[distList.index(min(distList))]])
            closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']

        #NOW WRITE THE ROW FOR THE ENHANCER TABLE
        newEnhancerLine = line[0:6]

        if byRefseq:
            newEnhancerLine.append(','.join(ROSE_utils.uniquify([x for x in overlappingGenes])))
            newEnhancerLine.append(','.join(ROSE_utils.uniquify([x for x in proximalGenes])))
            #print newEnhancerLine
            #print len(allEnhancerGenes)
            #print distList
            closestGene = allEnhancerGenes[distList.index(min(distList))]
            newEnhancerLine.append(closestGene)
        else:
            newEnhancerLine.append(','.join(ROSE_utils.uniquify([startDict[x]['name'] for x in overlappingGenes])))
            newEnhancerLine.append(','.join(ROSE_utils.uniquify([startDict[x]['name'] for x in proximalGenes])))
            closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']
            newEnhancerLine.append(closestGene)


        newEnhancerLine += line[-2:]
        enhancerToGeneTable.append(newEnhancerLine)
        #Now grab all overlapping and proximal genes for the gene ordered table

        overallGeneList +=overlappingGenes
        for refID in overlappingGenes:
            geneDict['overlapping'][refID].append(enhancerString)
            geneDict['enhancerString'][enhancerString].append(enhancerSignal)
        

        overallGeneList+=proximalGenes
        for refID in proximalGenes:
            geneDict['proximal'][refID].append(enhancerString)
            geneDict['enhancerString'][enhancerString].append(enhancerSignal)


    #End loop through
    
    #Make table by gene
    overallGeneList = ROSE_utils.uniquify(overallGeneList)
    

    nameOrder = ROSE_utils.order([startDict[x]['name'] for x in overallGeneList])

    usedNames = []

    for i in nameOrder:
        refID = overallGeneList[i]
        geneName = startDict[refID]['name']
        if usedNames.count(geneName) > 0 and uniqueGenes == True:

            continue
        else:
            usedNames.append(geneName)
        
        proxEnhancers = geneDict['proximal'][refID] + geneDict['overlapping'][refID]
        
        newLine = [geneName,refID,','.join(proxEnhancers)]

        
        for eachEnhancer in proxEnhancers:
            signalWithGenes.append([geneName,refID,eachEnhancer,ROSE_utils.uniquify(geneDict['enhancerString'][eachEnhancer])[0]])
        geneToEnhancerTable.append(newLine)

    #re-sort enhancerToGeneTable

    enhancerOrder = ROSE_utils.order([int(line[-2]) for line in enhancerToGeneTable[1:]])
    sortedTable = [enhancerToGeneTable[0]]
    for i in enhancerOrder:
        sortedTable.append(enhancerToGeneTable[(i+1)])

    return sortedTable,geneToEnhancerTable,signalWithGenes



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
                      help = "Enter the genome build (MM9,MM8,HG18,HG19,HG38)")
    parser.add_option("--custom", dest="custom_genome", default=None,
                      help = "Enter the custom genome annotation .ucsc")

    #optional flags
    parser.add_option("-l","--list", dest="geneList",nargs = 1, default=None,
                      help = "Enter a gene list to filter through")
    parser.add_option("-o","--out", dest="out",nargs = 1, default=None,
                      help = "Enter an output folder. Default will be same folder as input file")
    parser.add_option("-r","--refseq",dest="refseq",action = 'store_true', default=False,
                      help = "If flagged will write output by refseq ID and not common name")
    parser.add_option("-c","--control",dest="control",action = 'store_true', default=False,
                      help = "If flagged will subtract input from sample signal")

    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()


    if not options.input or not (options.genome or options.custom_genome):

        parser.print_help()
        exit()

    #GETTING THE INPUT
    enhancerFile = options.input

    #making the out folder if it doesn't exist
    if options.out:
        outFolder = ROSE_utils.formatFolder(options.out,True)
    else:
        outFolder = '/'.join(enhancerFile.split('/')[0:-1]) + '/'


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

    #GETTING THE GENOME
    if options.custom_genome:
        annotFile = options.custom_genome
        print('USING CUSTOM GENOME %s AS THE GENOME FILE' % options.custom_genome)
    else:
        genome = options.genome
        annotFile = genomeDict[genome.upper()]
        print('USING %s AS THE GENOME' % genome)

    #GETTING THE TRANSCRIBED LIST
    if options.geneList:

        transcribedFile = options.geneList
    else:
        transcribedFile = ''

    enhancerToGeneTable,geneToEnhancerTable,withGenesTable = mapEnhancerToGene(annotFile,enhancerFile, uniqueGenes=True, byRefseq=options.refseq, subtractInput=options.control, transcribedFile=transcribedFile)

    #Writing enhancer output
    enhancerFileName = enhancerFile.split('/')[-1].split('.')[0]

    #writing the enhancer table
    out1 = '%s%s_REGION_TO_GENE.txt' % (outFolder,enhancerFileName)
    ROSE_utils.unParseTable(enhancerToGeneTable,out1,'\t')

    #writing the gene table
    out2 = '%s%s_GENE_TO_REGION.txt' % (outFolder,enhancerFileName)
    ROSE_utils.unParseTable(geneToEnhancerTable,out2,'\t')

    out3 = '%s%s.table_withGENES.txt' % (outFolder,enhancerFileName)
    ROSE_utils.unParseTable(withGenesTable,out3,'\t')

if __name__ == "__main__":
    main()
