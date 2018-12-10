"""
THIS IS THE WORK DONE BY BY THE CHIMERA GROUP FOR THE BIOINF 575 FINAL

ASSIGNMENT DUE: 12/9/2018

-This program is meant to take in up to 4 .tsv files produced from STAR-Fusion and
 perform analysis on the chimeric fusion data produced from the STAR-Fusion programself.

-The program is only command line executable. In order to execute the program first go into
 whatever command line interface you use and change your directory to the one where this program
 is located. Second, use the following format  in order to run the program and produce the analysis
    $ python roazda_project.py (path to cell line 1's fusion data file) (path to cell line 2's fusion data file)
    Example:$ python roazda_project.py star-fusion.fusion_predictions_A172.tsv star-fusion.fusion_predictions_A375.tsv

-Input dependencies
    +in order for each table to reflect the cell line names properly, the title of the input
     files MUST BE in the format that STAR-Fusion produces them.
        IE: files must be in format of
                star-fusion.fusion_predictions_Cell-Line.tsv
    +FOR BRU-SEQ v BRU-CHASE data comparisons, IF ONLY 2 INPUT FILES ARE PROVIDED THE
     PROGRAM ASSUMES THAT YOU ARE DOING BRU-SEQ v BRU-CHASE ON THE SAME CELL LINE

REQUIRES:
    -at most 4 .tsv files that follow the format of the STAR-Fusion files found on
     the STAR-Fusion wiki (https://github.com/STAR-Fusion/STAR-Fusion/wiki)
    -providing only 2 .tsv files will have the program assume that you are attempting
     to do BRU-Chase v BRU-SEQ data anlysis
MODIFIES:
    -nothing
EFFECTS:
    -2 possible outputs from the file
        -STAR-DATA.txt: this file will include comparisons between 3-4 fusion files.
         Within each file is 4 tables
            +Cell Line Collection Summary
            +Frequent Genes and Their Relations
                -will display genes that appear AT LEAST 2X in at least 1 cell line
            +Fusions seen in 2 cell lines
            +Fusions seen in 3+ cell lines
        -STAR-Bru-Chase-Data.txt: this file will show comparisons between 2 files assuming that
         one file is BRU-Chase-Seq data and one is BRU-Seq data. Contains 3 tables
            +Cell Line Collection Summary
            +Frequent Genes and Their Relations
            +Fusions Shared
"""
import os
import sys
from statistics import stdev
from beautifultable import BeautifulTable
from collections import Counter


def main():
    #read in file names
    fileNameCellA = sys.argv[1]
    fileNameCellB = sys.argv[2]

    #collect names of the cell lines from each command line input
    #will contain a list of the filename when split at the _, name of cell line should be 3rd object (@ [2])
    nameCellLineA = sys.argv[1].strip().split('_')
    nameCellLineB = sys.argv[2].strip().split('_')


    #initialize all necessary containers for comparison for cell line A
    fusionNameCellA = set()
    leftGeneCellA = dict()
    leftGeneBreakCellA = dict()
    rightGeneCellA = dict()
    rightGeneBreakCellA = dict()
    annotsCellA = dict()
    leftAndRightGenesA = list()

    #all necessary containers for comparison for cell line B
    fusionNameCellB = set()
    leftGeneCellB = dict()
    leftGeneBreakCellB = dict()
    rightGeneCellB = dict()
    rightGeneBreakCellB = dict()
    annotsCellB = dict()
    leftAndRightGenesB = list()

    #all necessary containers for comparison for cell line C if it exists
    fusionNameCellC = set()
    leftGeneCellC = dict()
    leftGeneBreakCellC = dict()
    rightGeneCellC = dict()
    rightGeneBreakCellC = dict()
    annotsCellC = dict()
    leftAndRightGenesC = list()

    #all necessary containers for comparison for cell line D if it exists
    fusionNameCellD = set()
    leftGeneCellD = dict()
    leftGeneBreakCellD = dict()
    rightGeneCellD = dict()
    rightGeneBreakCellD = dict()
    annotsCellD = dict()
    leftAndRightGenesD = list()

    #create a list that will contain all genes from all cell lines
    totalGeneList = list()

    #read thru each file and append/extend/update each container as necessary
    if os.path.exists(fileNameCellA):
        with open(fileNameCellA, 'r') as infile:
            for line in infile:
                #should ignore the immediate header row
                if line.startswith("#"): continue
                #create a list of each string in a row
                parts = line.strip().split("\t")
                #split fusion entry to get left and right gene
                left_rightGene = parts[0].strip().split('--')
                left_break = parts[5].strip().split(':')
                right_break = parts[7].strip().split(':')
                #add each element into its respective container
                fusionNameCellA.add(parts[0])
                leftGeneCellA.update({parts[0]: left_rightGene[0]})
                leftGeneBreakCellA.update({parts[0]: left_break[1]})
                rightGeneCellA.update({parts[0]: left_rightGene[1]})
                rightGeneBreakCellA.update({parts[0]: right_break[1]})
                annotsCellA.update({parts[0]: parts[16]})
                leftAndRightGenesA.extend((left_rightGene[0], left_rightGene[1]))


    if os.path.exists(fileNameCellB):
        with open(fileNameCellB, 'r') as infile:
            for line in infile:
                #should ignore the immediate header row
                if line.startswith("#"): continue
                #create a list of each string in a row
                parts = line.strip().split("\t")
                #split fusion entry to get left and right gene
                left_rightGene = parts[0].split('--')
                left_break = parts[5].strip().split(':')
                right_break = parts[7].strip().split(':')
                #add each element into its respective container
                fusionNameCellB.add(parts[0])
                leftGeneCellB.update({parts[0]: left_rightGene[0]})
                leftGeneBreakCellB.update({parts[0]: left_break[1]})
                rightGeneCellB.update({parts[0]: left_rightGene[1]})
                rightGeneBreakCellB.update({parts[0]: right_break[1]})
                annotsCellB.update({parts[0]: parts[16]})
                leftAndRightGenesB.extend((left_rightGene[0], left_rightGene[1]))
    #if we are given 4 cell lines, then we know we arent handling bru-chase/seq and we should do proper analysis
    if len(sys.argv) > 3:
        fileNameCellC = sys.argv[3]
        fileNameCellD = sys.argv[4]
        nameCellLineC = sys.argv[3].strip().split('_')
        nameCellLineD = sys.argv[4].strip().split('_')
        if os.path.exists(fileNameCellC):
            with open(fileNameCellC, 'r') as infile:
                for line in infile:
                    #should ignore the immediate header row
                    if line.startswith("#"): continue
                    #create a list of each string in a row
                    parts = line.strip().split("\t")
                    #split fusion entry to get left and right gene
                    left_rightGene = parts[0].split('--')
                    left_break = parts[5].strip().split(':')
                    right_break = parts[7].strip().split(':')
                    #add each element into its respective container
                    fusionNameCellC.add(parts[0])
                    leftGeneCellC.update({parts[0]: left_rightGene[0]})
                    leftGeneBreakCellC.update({parts[0]: left_break[1]})
                    rightGeneCellC.update({parts[0]: left_rightGene[1]})
                    rightGeneBreakCellC.update({parts[0]: right_break[1]})
                    annotsCellC.update({parts[0]: parts[16]})
                    leftAndRightGenesC.extend((left_rightGene[0], left_rightGene[1]))
        if os.path.exists(fileNameCellD):
            with open(fileNameCellD, 'r') as infile:
                for line in infile:
                    #should ignore the immediate header row
                    if line.startswith("#"): continue
                    #create a list of each string in a row
                    parts = line.strip().split("\t")
                    #split fusion entry to get left and right gene
                    left_rightGene = parts[0].split('--')
                    left_break = parts[5].strip().split(':')
                    right_break = parts[7].strip().split(':')
                    #add each element into its respective container
                    fusionNameCellD.add(parts[0])
                    leftGeneCellD.update({parts[0]: left_rightGene[0]})
                    leftGeneBreakCellD.update({parts[0]: left_break[1]})
                    rightGeneCellD.update({parts[0]: left_rightGene[1]})
                    rightGeneBreakCellD.update({parts[0]: right_break[1]})
                    annotsCellD.update({parts[0]: parts[16]})
                    leftAndRightGenesD.extend((left_rightGene[0], left_rightGene[1]))

        #add to each set the fusions that are shared among the cell lines being compared
            #done using the set.intersection function to isolate those shared fusions
        totalGeneList = leftAndRightGenesA + leftAndRightGenesB + leftAndRightGenesC + leftAndRightGenesD
        AtoBNameConnections = set.intersection(fusionNameCellA, fusionNameCellB)
        AtoCNameConnections = set.intersection(fusionNameCellA, fusionNameCellC)
        AtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellD)
        BtoCNameConnections = set.intersection(fusionNameCellB, fusionNameCellC)
        BtoDNameConnections = set.intersection(fusionNameCellB, fusionNameCellD)
        CtoDNameConnections = set.intersection(fusionNameCellC, fusionNameCellD)
        AtoBtoCNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC)
        AtoBtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellD)
        AtoCtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellC, fusionNameCellD)
        BtoCtoDNameConnections = set.intersection(fusionNameCellB, fusionNameCellC, fusionNameCellD)
        AtoBtoCtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC, fusionNameCellD)

        #create a Counter that will organize the genes into most frequent
            #-the counter is a dictionary where each key is each gene and the value will be its frequency
        mostComGenesA = Counter(leftAndRightGenesA)
        mostComGenesB = Counter(leftAndRightGenesB)
        mostComGenesC = Counter(leftAndRightGenesC)
        mostComGenesD = Counter(leftAndRightGenesD)

        #make sets of the intersection between different cell lines left and right genes list to find the common genes between cell lines
        comGenesA_B = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesB))
        comGenesA_C = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesC))
        comGenesA_D = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesD))
        comGenesB_C = set.intersection(set(leftAndRightGenesB), set(leftAndRightGenesC))
        comGenesB_D = set.intersection(set(leftAndRightGenesB), set(leftAndRightGenesD))
        comGenesC_D = set.intersection(set(leftAndRightGenesC), set(leftAndRightGenesD))
        comGenesA_B_C = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesB), set(leftAndRightGenesC))
        comGenesA_B_D =set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesB), set(leftAndRightGenesD))
        comGenesA_C_D = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesC), set(leftAndRightGenesD))
        comGenesB_C_D = set.intersection(set(leftAndRightGenesB), set(leftAndRightGenesC), set(leftAndRightGenesD))
        comGenesA_B_C_D = set.intersection(set(leftAndRightGenesB), set(leftAndRightGenesA), set(leftAndRightGenesC), set(leftAndRightGenesD))


        #create table to show number of fusions shared and number of genes shared between cell lines
            #beautifultable is a python extension that allows our data to be printed to the terminal in a lovely clear fashion
        FusionNumberTable = BeautifulTable(max_width = 200)
        FusionNumberTable.column_headers = ["Connection", "# of Fusions Shared", "# of Genes Shared"]
        FusionNumberTable.append_row(["A-B", len(AtoBNameConnections), len(comGenesA_B)])
        FusionNumberTable.append_row(["A-C", len(AtoCNameConnections), len(comGenesA_C)])
        FusionNumberTable.append_row(["A-D", len(AtoDNameConnections), len(comGenesA_D)])
        FusionNumberTable.append_row(["B-C", len(BtoCNameConnections), len(comGenesB_C)])
        FusionNumberTable.append_row(["B-D", len(BtoDNameConnections), len(comGenesB_D)])
        FusionNumberTable.append_row(["C-D", len(CtoDNameConnections), len(comGenesC_D)])
        FusionNumberTable.append_row(["A-B-C", len(AtoBtoCNameConnections), len(comGenesA_B_C)])
        FusionNumberTable.append_row(["A-B-D", len(AtoBtoDNameConnections), len(comGenesA_B_D)])
        FusionNumberTable.append_row(["A-C-D", len(AtoCtoDNameConnections), len(comGenesA_C_D)])
        FusionNumberTable.append_row(["B-C-D", len(BtoCtoDNameConnections), len(comGenesB_C_D)])
        FusionNumberTable.append_row(["A-B-C-D", len(AtoBtoCtoDNameConnections), len(comGenesA_B_C_D)])
        FusionNumberTable.left_border_char = '|'
        FusionNumberTable.right_border_char = '|'
        FusionNumberTable.top_border_char = '~'
        FusionNumberTable.bottom_border_char = '~'
        FusionNumberTable.header_separator_char = '='
        FusionNumberTable.row_separator_char = ''
        FusionNumberTable.intersection_char = ''
        FusionNumberTable.column_separator_char = ':'

        #Create the tables for each gene fusion connection
        FusionTable = BeautifulTable(max_width = 200)
        FusionTable.column_headers = ["Connection","Fusion Name", "Left Gene", "Right Gene", "Difference in Left Gene Breakpoint", "Difference in Right Gene Breakpoint", "Annotation of Fusion"]
        #get each left and right gene from each respective fusion dictionary, subtract the difference between the breakpoints
        for fusion in AtoBNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellB[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellB[fusion])
            FusionTable.append_row(["A-B", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion]])
        for fusion in AtoCNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellC[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellC[fusion])
            FusionTable.append_row(["A-C", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion]])
        for fusion in AtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row(["A-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion]])
        for fusion in BtoCNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakDiff = int(leftGeneBreakCellB[fusion]) - int(leftGeneBreakCellC[fusion])
            rightBreakDiff = int(rightGeneBreakCellB[fusion]) - int(rightGeneBreakCellC[fusion])
            FusionTable.append_row(["B-C", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellB[fusion]])
        for fusion in BtoDNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakDiff = int(leftGeneBreakCellB[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellB[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row(["B-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellB[fusion]])
        for fusion in CtoDNameConnections:
            leftGene = leftGeneCellC[fusion]
            rightGene = rightGeneCellC[fusion]
            leftBreakDiff = int(leftGeneBreakCellC[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellC[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row(["C-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellC[fusion]])
        #format the table to make it look all pretty :D
        FusionTable.left_border_char = '|'
        FusionTable.right_border_char = '|'
        FusionTable.top_border_char = '~'
        FusionTable.bottom_border_char = '~'
        FusionTable.header_separator_char = '='
        FusionTable.row_separator_char = ''
        FusionTable.intersection_char = ''
        FusionTable.column_separator_char = ':'
        #create table to describe connections of fusions in 3+ cell lines
        FusionTableThreePlus = BeautifulTable(max_width = 200)
        FusionTableThreePlus.column_headers = ["Fusion", "Connection", "Left Gene", "Right Gene", "Left Break St Dev", "Right Break St Dev", "Annotation"]
        #similar procedure to the previous table but instead of finding the difference, we use the standard deviation function
        for fusion in AtoBtoCNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [int(leftGeneBreakCellA[fusion]), int(leftGeneBreakCellB[fusion]), int(leftGeneBreakCellC[fusion])]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [int(rightGeneBreakCellA[fusion]), int(rightGeneBreakCellB[fusion]), int(rightGeneBreakCellC[fusion])]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row([fusion, 'A-B-C', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion]])
        for fusion in AtoBtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [int(leftGeneBreakCellA[fusion]), int(leftGeneBreakCellB[fusion]), int(leftGeneBreakCellD[fusion])]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [int(rightGeneBreakCellA[fusion]), int(rightGeneBreakCellB[fusion]), int(rightGeneBreakCellD[fusion])]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row([fusion, 'A-B-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion]])
        for fusion in AtoCtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [int(leftGeneBreakCellA[fusion]), int(leftGeneBreakCellC[fusion]), int(leftGeneBreakCellD[fusion])]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [int(rightGeneBreakCellA[fusion]), int(rightGeneBreakCellC[fusion]), int(rightGeneBreakCellD[fusion])]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row([fusion, 'A-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion]])
        for fusion in BtoCtoDNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakSample = [int(leftGeneBreakCellB[fusion]), int(leftGeneBreakCellC[fusion]), int(leftGeneBreakCellD[fusion])]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [int(rightGeneBreakCellB[fusion]), int(rightGeneBreakCellC[fusion]), int(rightGeneBreakCellD[fusion])]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row([fusion, 'B-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion]])
        for fusion in AtoBtoCtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [int(leftGeneBreakCellA[fusion]), int(leftGeneBreakCellB[fusion]), int(leftGeneBreakCellC[fusion]), int(leftGeneBreakCellD[fusion])]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [int(rightGeneBreakCellA[fusion]), int(rightGeneBreakCellB[fusion]), int(rightGeneBreakCellC[fusion]), int(rightGeneBreakCellD[fusion])]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row([fusion, 'A-B-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion]])
        FusionTableThreePlus.left_border_char = '|'
        FusionTableThreePlus.right_border_char = '|'
        FusionTableThreePlus.top_border_char = '~'
        FusionTableThreePlus.bottom_border_char = '~'
        FusionTableThreePlus.header_separator_char = '='
        FusionTableThreePlus.row_separator_char = ''
        FusionTableThreePlus.intersection_char = ''
        FusionTableThreePlus.column_separator_char = ':'
        #create table to describe the genes seen multiple times within cell lines
        FusionTableGeneCommonalities = BeautifulTable(max_width = 200)
        FusionTableGeneCommonalities.column_headers = ["Common Gene", "Occurrences in Cell Line A", "Occurrences in Cell Line B", "Occurrences in Cell Line C", "Occurrences in Cell Line D"]
        #if a gene occurs more than twice in one cell line then it will be added to the table
        for gene in set(totalGeneList):
            a_occur = 0
            b_occur = 0
            c_occur = 0
            d_occur = 0
            if gene in mostComGenesA:
                a_occur = mostComGenesA[gene]
            if gene in mostComGenesB:
                b_occur = mostComGenesB[gene]
            if gene in mostComGenesC:
                c_occur = mostComGenesC[gene]
            if gene in mostComGenesD:
                d_occur = mostComGenesD[gene]
            if a_occur > 1 or b_occur > 1 or c_occur > 1 or d_occur >1:
                FusionTableGeneCommonalities.append_row([gene, a_occur, b_occur, c_occur, d_occur])
        FusionTableGeneCommonalities.left_border_char = '|'
        FusionTableGeneCommonalities.right_border_char = '|'
        FusionTableGeneCommonalities.top_border_char = '~'
        FusionTableGeneCommonalities.bottom_border_char = '~'
        FusionTableGeneCommonalities.header_separator_char = '='
        FusionTableGeneCommonalities.row_separator_char = ''
        FusionTableGeneCommonalities.intersection_char = ''
        FusionTableGeneCommonalities.column_separator_char = ':'

        #some formating of the table before we print
        FusionTable.width_exceed_policy = BeautifulTable.WEP_WRAP
        FusionTableThreePlus.width_exceed_policy = BeautifulTable.WEP_WRAP
        FusionTableGeneCommonalities.width_exceed_policy = BeautifulTable.WEP_WRAP
        FusionNumberTable.width_exceed_policy = BeautifulTable.WEP_WRAP
        FusionTable.sort("Fusion Name")
        FusionTableThreePlus.sort("Fusion")
        FusionTableGeneCommonalities.sort("Common Gene")
        FusionNumberTable.sort("# of Fusions Shared", reverse = 1)
        uniqueGenesA = set(list(leftGeneCellA.values()) + list(rightGeneCellA.values()))
        uniqueGenesB = set(list(leftGeneCellB.values()) + list(rightGeneCellB.values()))
        uniqueGenesC = set(list(leftGeneCellC.values()) + list(rightGeneCellD.values()))
        uniqueGenesD = set(list(leftGeneCellD.values()) + list(rightGeneCellD.values()))



        file = open("Star-Data.txt", "w")
        #start printing everything
        file.write("A refers to %s, B refers to %s, C refers to %s, D refers to %s" % (nameCellLineA[2], nameCellLineB[2], nameCellLineC[2], nameCellLineD[2]))
        file.write('\n')
        file.write('\n')
        file.write("Number of unique genes in {} = {} \n".format(nameCellLineA[2], len(uniqueGenesA)))
        file.write("Number of unique genes in {} = {} \n".format(nameCellLineB[2], len(uniqueGenesB)))
        file.write("Number of unique genes in {} = {} \n".format(nameCellLineC[2], len(uniqueGenesC)))
        file.write("Number of unique genes in {} = {} \n".format(nameCellLineD[2], len(uniqueGenesD)))
        file.write('\n')
        file.write("Cell Line Collection Summary")
        file.write('\n')
        file.write(str(FusionNumberTable))

        for i in range(3):
            file.write('\n')
        file.write("Frequent Genes and Their Relations")
        file.write('\n')
        file.write(str(FusionTableGeneCommonalities))

        for i in range(3):
            file.write('\n')
        file.write("Fusions Seen in 2 Cell Lines")
        file.write('\n')
        file.write(str(FusionTable))

        for i in range(3):
            file.write('\n')
        file.write("Fusions Seen in 3+ Cell Lines")
        file.write('\n')
        file.write(str(FusionTableThreePlus))

        file.close()

    #dealing with the bru-chase/seq data when only two folders are given on the command line
    else:
        totalGeneList = leftAndRightGenesA + leftAndRightGenesB
        AtoBNameConnections = set.intersection(fusionNameCellA, fusionNameCellB)
        mostComGenesA = Counter(leftAndRightGenesA)
        mostComGenesB = Counter(leftAndRightGenesB)
        comGenesA_B = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesB))
        #make its own general conection table
        FusionNumberTable = BeautifulTable(max_width = 200)
        FusionNumberTable.column_headers = ["Connection", "# of Fusions Shared", "# of Genes Shared"]
        FusionNumberTable.append_row(["A-B", len(AtoBNameConnections), len(comGenesA_B)])
        FusionNumberTable.left_border_char = '|'
        FusionNumberTable.right_border_char = '|'
        FusionNumberTable.top_border_char = '~'
        FusionNumberTable.bottom_border_char = '~'
        FusionNumberTable.header_separator_char = '='
        FusionNumberTable.row_separator_char = ''
        FusionNumberTable.intersection_char = ''
        FusionNumberTable.column_separator_char = ':'

        #make its own connection tables
        FusionTable = BeautifulTable(max_width = 200)
        FusionTable.column_headers = ["Fusion Name", "Left Gene", "Right Gene", "Difference in Left Gene Breakpoint", "Difference in Right Gene Breakpoint", "Annotation of Fusion"]
        for fusion in AtoBNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellB[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellB[fusion])
            FusionTable.append_row([fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion]])
        FusionTable.left_border_char = '|'
        FusionTable.right_border_char = '|'
        FusionTable.top_border_char = '~'
        FusionTable.bottom_border_char = '~'
        FusionTable.header_separator_char = '='
        FusionTable.row_separator_char = ''
        FusionTable.intersection_char = ''
        FusionTable.column_separator_char = ':'
        #make a table showing commonalities in gene
        FusionTableGeneCommonalities = BeautifulTable(max_width = 200)
        FusionTableGeneCommonalities.column_headers = ["Common Gene", "Occurrences in Cell Line A", "Occurrences in Cell Line B"]
        for gene in set(totalGeneList):
            a_occur = 0
            b_occur = 0
            if gene in mostComGenesA:
                a_occur = mostComGenesA[gene]
            if gene in mostComGenesB:
                b_occur = mostComGenesB[gene]
            if a_occur > 1 or b_occur > 1:
                FusionTableGeneCommonalities.append_row([gene, a_occur, b_occur])
        FusionTableGeneCommonalities.left_border_char = '|'
        FusionTableGeneCommonalities.right_border_char = '|'
        FusionTableGeneCommonalities.top_border_char = '~'
        FusionTableGeneCommonalities.bottom_border_char = '~'
        FusionTableGeneCommonalities.header_separator_char = '='
        FusionTableGeneCommonalities.row_separator_char = ''
        FusionTableGeneCommonalities.intersection_char = ''
        FusionTableGeneCommonalities.column_separator_char = ':'

        file = open("Star-Bru-Chase-Data.txt", "w")

        file.write("A refers to %s, B refers to %s of %s cell line" % (nameCellLineA[3], nameCellLineB[3], nameCellLineA[2]))
        for i in range(3):
            file.write('\n')
        file.write("Cell Line Collection Summary")
        file.write('\n')
        file.write(str(FusionNumberTable))

        for i in range(3):
            file.write('\n')
        file.write("Frequent Genes and Their Relations")
        file.write('\n')
        file.write(str(FusionTableGeneCommonalities))

        for i in range(3):
            file.write('\n')
        file.write("Fusions Seen in 2 Cell Lines")
        file.write('\n')
        file.write(str(FusionTable))

        file.close()
if __name__ == '__main__':
    main()
