import os
import sys
from statistics import stdev
from beautifultable import BeautifulTable
from collections import Counter


def main():
    #read in all file names
    fileNameCellA = sys.argv[1]
    fileNameCellB = sys.argv[2]
    fileNameCellC = sys.argv[3]
    fileNameCellD = sys.argv[4]


    #all necessary containers for comparison for cell line A
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

    totalGeneList = list()

    #read thru each file and append/extend/update each container as necessary
    if os.path.exists(fileNameCellA):
        with open(fileNameCellA, 'r') as infile:
            for line in infile:
                #should ignore the immediate header files
                if line.startswith("#"): continue
                #create a list of each string in a row
                parts = line.strip().split("\t")
                #split fusion entry to get left and right gene
                left_rightGene = parts[0].strip().split('--')
                left_break = parts[5].strip().split(':')
                right_break = parts[7].strip().split(':')
                #add each element into its respective dictionary
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
                #should ignore the immediate header files
                if line.startswith("#"): continue
                #create a list of each string in a row
                parts = line.strip().split("\t")
                #split fusion entry to get left and right gene
                left_rightGene = parts[0].split('--')
                left_break = parts[5].strip().split(':')
                right_break = parts[7].strip().split(':')
                #add each element into its respective dictionary
                fusionNameCellB.add(parts[0])
                leftGeneCellB.update({parts[0]: left_rightGene[0]})
                leftGeneBreakCellB.update({parts[0]: left_break[1]})
                rightGeneCellB.update({parts[0]: left_rightGene[1]})
                rightGeneBreakCellB.update({parts[0]: right_break[1]})
                annotsCellB.update({parts[0]: parts[16]})
                leftAndRightGenesB.extend((left_rightGene[0], left_rightGene[1]))

    if sys.argv[3] and sys.argv[4]:
        if os.path.exists(fileNameCellC):
            with open(fileNameCellC, 'r') as infile:
                for line in infile:
                    #should ignore the immediate header files
                    if line.startswith("#"): continue
                    #create a list of each string in a row
                    parts = line.strip().split("\t")
                    #split fusion entry to get left and right gene
                    left_rightGene = parts[0].split('--')
                    left_break = parts[5].strip().split(':')
                    right_break = parts[7].strip().split(':')
                    #add each element into its respective dictionary
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
                    #should ignore the immediate header files
                    if line.startswith("#"): continue
                    #create a list of each string in a row
                    parts = line.strip().split("\t")
                    #split fusion entry to get left and right gene
                    left_rightGene = parts[0].split('--')
                    left_break = parts[5].strip().split(':')
                    right_break = parts[7].strip().split(':')
                    #add each element into its respective dictionary
                    fusionNameCellD.add(parts[0])
                    leftGeneCellD.update({parts[0]: left_rightGene[0]})
                    leftGeneBreakCellD.update({parts[0]: left_break[1]})
                    rightGeneCellD.update({parts[0]: left_rightGene[1]})
                    rightGeneBreakCellD.update({parts[0]: right_break[1]})
                    annotsCellD.update({parts[0]: parts[16]})
                    leftAndRightGenesD.extend((left_rightGene[0], left_rightGene[1]))

    #add to each set the fusions that are shared among the cell lines being compared
        totalGeneList = leftAndRightGenesA + leftAndRightGenesB + leftAndRightGenesC + leftAndRightGenesD
        AtoBNameConnections = set.intersection(fusionNameCellA, fusionNameCellB)
        AtoCNameConnections = set.intersection(fusionNameCellA, fusionNameCellC)
        AtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellD)
        BtoCNameConnections = set.intersection(fusionNameCellB, fusionNameCellC)
        BtoDNameConnections = set.intersection(fusionNameCellB, fusionNameCellD)
        CtoDNameConnections = set.intersection(fusionNameCellD, fusionNameCellD)
        AtoBtoCNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC)
        AtoBtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellD)
        AtoCtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellC, fusionNameCellD)
        BtoCtoDNameConnections = set.intersection(fusionNameCellB, fusionNameCellC, fusionNameCellD)
        AtoBtoCtoDNameConnections = set.intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC, fusionNameCellD)

    #create a Counter that will organize the genes into most frequent
    #the counter is a dictionary where each key is each gene and the value will be its frequency
        mostComGenesA = Counter(leftAndRightGenesA)
        mostComGenesB = Counter(leftAndRightGenesB)
        mostComGenesC = Counter(leftAndRightGenesC)
        mostComGenesD = Counter(leftAndRightGenesD)
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

    #create table to show number of fusions shared
        FusionNumberTable = BeautifulTable()
        FusionNumberTable.column_headers("Connection", "# of Fusions Shared", "# of Genes Shared")
        FusionNumberTable.append_row("A-B", len(AtoBNameConnections), len(comGenesA_B))
        FusionNumberTable.append_row("A-C", len(AtoCNameConnections), len(comGenesA_C))
        FusionNumberTable.append_row("A-D", len(AtoDNameConnections), len(comGenesA_D))
        FusionNumberTable.append_row("B-C", len(BtoCNameConnections), len(comGenesB_C))
        FusionNumberTable.append_row("B-D", len(BtoDNameConnections), len(comGenesB_D))
        FusionNumberTable.append_row("C-D", len(CtoDNameConnections), len(comGenesC_D))
        FusionNumberTable.append_row("A-B-C", len(AtoBtoCNameConnections), len(comGenesA_B_C))
        FusionNumberTable.append_row("A-B-D", len(AtoBtoDNameConnections), len(comGenesA_B_D))
        FusionNumberTable.append_row("A-C-D", len(AtoCtoDNameConnections), len(comGenesA_C_D))
        FusionNumberTable.append_row("B-C-D", len(BtoCtoDNameConnections), len(comGenesB_C_D))
        FusionNumberTable.append_row("A-B-C-D", len(AtoBtoCtoDNameConnections), len(comGenesA_B_C_D))
        FusionNumberTable.left_border_char = 'l'
        FusionNumberTable.right_border_char = 'l'
        FusionNumberTable.top_border_char = '~'
        FusionNumberTable.bottom_border_char = '~'
        FusionNumberTable.header_separator_char = '='
        FusionNumberTable.row_separator_char = ''
        FusionNumberTable.intersection_char = ''
        FusionNumberTable.column_separator_char = ':'


    #Create the tables for each gene fusion connection
    #beautifultable is a python extension that allows our data to be printed to the terminal in a lovely clear fashion
        FusionTable = BeautifulTable()
        FusionTable.column_headers("Connection","Fusion Name", "Left Gene", "Right Gene", "Difference in Left Gene Breakpoint", "Difference in Right Gene Breakpoint", "Annotation of Fusion")
        for fusion in AtoBNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellB[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellB[fusion])
            FusionTable.append_row("A-B", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion])
        for fusion in AtoCNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellC[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellC[fusion])
            FusionTable.append_row("A-C", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion])
        for fusion in AtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row("A-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion])
        for fusion in BtoCNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakDiff = int(leftGeneBreakCellB[fusion]) - int(leftGeneBreakCellC[fusion])
            rightBreakDiff = int(rightGeneBreakCellB[fusion]) - int(rightGeneBreakCellC[fusion])
            FusionTable.append_row("B-C", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellB[fusion])
        for fusion in BtoDNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakDiff = int(leftGeneBreakCellB[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellB[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row("B-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellB[fusion])
        for fusion in CtoDNameConnections:
            leftGene = leftGeneCellC[fusion]
            rightGene = rightGeneCellC[fusion]
            leftBreakDiff = int(leftGeneBreakCellC[fusion]) - int(leftGeneBreakCellD[fusion])
            rightBreakDiff = int(rightGeneBreakCellC[fusion]) - int(rightGeneBreakCellD[fusion])
            FusionTable.append_row("C-D", fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellC[fusion])
        FusionTable.left_border_char = 'l'
        FusionTable.right_border_char = 'l'
        FusionTable.top_border_char = '~'
        FusionTable.bottom_border_char = '~'
        FusionTable.header_separator_char = '='
        FusionTable.row_separator_char = ''
        FusionTable.intersection_char = ''
        FusionTable.column_separator_char = ':'

        #create table to describe connections of fusions in 3+ cell lines
        FusionTableThreePlus = BeautifulTable()
        FusionTableThreePlus.column_headers("Fusion", "Connection", "Left Gene", "Right Gene", "Left Break St Dev", "Right Break St Dev", "Annotation")
        for fusion in AtoBtoCNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [leftGeneBreakCellA[fusion], leftGeneBreakCellB[fusion], leftGeneBreakCellC[fusion]]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [rightGeneBreakCellA[fusion], rightGeneBreakCellB[fusion], rightGeneBreakCellC[fusion]]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row(fusion, 'A-B-C', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion])
        for fusion in AtoBtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [leftGeneBreakCellA[fusion], leftGeneBreakCellB[fusion], leftGeneBreakCellD[fusion]]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [rightGeneBreakCellA[fusion], rightGeneBreakCellB[fusion], rightGeneBreakCellD[fusion]]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row(fusion, 'A-B-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion])
        for fusion in AtoCtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [leftGeneBreakCellA[fusion], leftGeneBreakCellC[fusion], leftGeneBreakCellD[fusion]]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [rightGeneBreakCellA[fusion], rightGeneBreakCellC[fusion], rightGeneBreakCellD[fusion]]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row(fusion, 'A-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion])
        for fusion in BtoCtoDNameConnections:
            leftGene = leftGeneCellB[fusion]
            rightGene = rightGeneCellB[fusion]
            leftBreakSample = [leftGeneBreakCellB[fusion], leftGeneBreakCellC[fusion], leftGeneBreakCellD[fusion]]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [rightGeneBreakCellB[fusion], rightGeneBreakCellC[fusion], rightGeneBreakCellD[fusion]]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row(fusion, 'B-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion])
        for fusion in AtoBtoCtoDNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakSample = [leftGeneBreakCellA[fusion], leftGeneBreakCellC[fusion], leftGeneBreakCellD[fusion], leftGeneBreakCellD]
            leftStDev = stdev(leftBreakSample)
            rightBreakSample = [rightGeneBreakCellA[fusion], rightGeneBreakCellC[fusion], rightGeneBreakCellD[fusion], rightGeneBreakCellD]
            rightStDev = stdev(rightBreakSample)
            FusionTableThreePlus.append_row(fusion, 'A-B-C-D', leftGene, rightGene, leftStDev, rightStDev, annotsCellA[fusion])
        FusionTableThreePlus.left_border_char = 'l'
        FusionTableThreePlus.right_border_char = 'l'
        FusionTableThreePlus.top_border_char = '~'
        FusionTableThreePlus.bottom_border_char = '~'
        FusionTableThreePlus.header_separator_char = '='
        FusionTableThreePlus.row_separator_char = ''
        FusionTableThreePlus.intersection_char = ''
        FusionTableThreePlus.column_separator_char = ':'


        #create table to describe the genes seen multiple times within cell lines
        FusionTableGeneCommonalities = BeautifulTable()
        FusionTableGeneCommonalities.column_headers("Common Gene", "Occurences in Cell Line A", "Occurences in Cell Line B", "Occurences in Cell Line C", "Occurences in Cell Line D")
        for gene in totalGeneList:
            a_occur = 0
            b_occur = 0
            c_occur = 0
            d_occur = 0
            if gene in mostComGenesA and mostComGenesA[gene] > 1:
                a_occur = mostComGenesA[gene]
            if gene in mostComGenesB and mostComGenesB[gene] > 1:
                b_occur = mostComGenesB[gene]
            if gene in mostComGenesC and mostComGenesC[gene] > 1:
                c_occur = mostComGenesC[gene]
            if gene in mostComGenesD and mostComGenesD[gene] > 1:
                d_occur = mostComGenesD[gene]
            if a_occur > 1 or b_occur > 1 or c_occur > 1 or d_occur >1:
                FusionTableGeneCommonalities.append_row(gene, a_occur, b_occur, c_occur, d_occur)
        FusionTableGeneCommonalities.left_border_char = 'l'
        FusionTableGeneCommonalities.right_border_char = 'l'
        FusionTableGeneCommonalities.top_border_char = '~'
        FusionTableGeneCommonalities.bottom_border_char = '~'
        FusionTableGeneCommonalities.header_separator_char = '='
        FusionTableGeneCommonalities.row_separator_char = ''
        FusionTableGeneCommonalities.intersection_char = ''
        FusionTableGeneCommonalities.column_separator_char = ':'

    #dealing with the bru-chase/seq data when only two folders are given on the command line
    else:
        totalGeneList = leftAndRightGenesA + leftAndRightGenesB
        AtoBNameConnections = set.intersection(fusionNameCellA, fusionNameCellB)
        mostComGenesA = Counter(leftAndRightGenesA)
        mostComGenesB = Counter(leftAndRightGenesB)
        comGenesA_B = set.intersection(set(leftAndRightGenesA), set(leftAndRightGenesB))
        #make its own general conection table
        FusionNumberTable = BeautifulTable()
        FusionNumberTable.column_headers("Connection", "# of Fusions Shared", "# of Genes Shared")
        FusionNumberTable.append_row("A-B", len(AtoBNameConnections), len(comGenesA_B))
        FusionNumberTable.left_border_char = 'l'
        FusionNumberTable.right_border_char = 'l'
        FusionNumberTable.top_border_char = '~'
        FusionNumberTable.bottom_border_char = '~'
        FusionNumberTable.header_separator_char = '='
        FusionNumberTable.row_separator_char = ''
        FusionNumberTable.intersection_char = ''
        FusionNumberTable.column_separator_char = ':'

        #make its own connection tables
        FusionTable = BeautifulTable()
        FusionTable.column_headers("Fusion Name", "Left Gene", "Right Gene", "Difference in Left Gene Breakpoint", "Difference in Right Gene Breakpoint", "Annotation of Fusion")
        for fusion in AtoBNameConnections:
            leftGene = leftGeneCellA[fusion]
            rightGene = rightGeneCellA[fusion]
            leftBreakDiff = int(leftGeneBreakCellA[fusion]) - int(leftGeneBreakCellB[fusion])
            rightBreakDiff = int(rightGeneBreakCellA[fusion]) - int(rightGeneBreakCellB[fusion])
            FusionTable.append_row(fusion, leftGene, rightGene, leftBreakDiff, rightBreakDiff, annotsCellA[fusion])
        FusionTable.left_border_char = 'l'
        FusionTable.right_border_char = 'l'
        FusionTable.top_border_char = '~'
        FusionTable.bottom_border_char = '~'
        FusionTable.header_separator_char = '='
        FusionTable.row_separator_char = ''
        FusionTable.intersection_char = ''
        FusionTable.column_separator_char = ':'
        #make a table showing commonalities in gene
        FusionTableGeneCommonalities = BeautifulTable()
        FusionTableGeneCommonalities.column_headers("Common Gene", "Occurences in Cell Line A", "Occurences in Cell Line B")
        for gene in totalGeneList:
            a_occur = 0
            b_occur = 0
            if gene in mostComGenesA and mostComGenesA[gene] > 1:
                a_occur = mostComGenesA[gene]
            if gene in mostComGenesB and mostComGenesB[gene] > 1:
                b_occur = mostComGenesB[gene]
            if a_occur > 1 or b_occur > 1 or c_occur > 1 or d_occur >1:
                FusionTableGeneCommonalities.append_row(gene, a_occur, b_occur)
        FusionTableGeneCommonalities.left_border_char = 'l'
        FusionTableGeneCommonalities.right_border_char = 'l'
        FusionTableGeneCommonalities.top_border_char = '~'
        FusionTableGeneCommonalities.bottom_border_char = '~'
        FusionTableGeneCommonalities.header_separator_char = '='
        FusionTableGeneCommonalities.row_separator_char = ''
        FusionTableGeneCommonalities.intersection_char = ''
        FusionTableGeneCommonalities.column_separator_char = ':'










if __name__ == '__main__':
    main()
