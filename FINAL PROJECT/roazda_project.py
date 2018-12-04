import os
import sys
from beautifultable import BeautifulTable
from collections import Counter


def main():
    fileNameCellA = sys.argv[1]
    fileNameCellB = sys.argv[2]
    fileNameCellC = sys.argv[3]
    fileNameCellD = sys.argv[4]


    #all necessary containers for comparison for cell line A
    fusionNameCellA = set()
    junctionReadCountCellA = dict()
    spanningFragCountCellA= dict()
    spliceTypeCellA = dict()
    leftGeneCellA = dict()
    leftGeneBreakCellA = dict()
    rightGeneCellA = dict()
    rightGeneBreakCellA = dict()
    annotsCellA = dict()
    leftAndRightGenesA = list()

    #all necessary containers for comparison for cell line B
    fusionNameCellB = set()
    junctionReadCountCellB = dict()
    spanningFragCountCellB= dict()
    spliceTypeCellB = dict()
    leftGeneCellB = dict()
    leftGeneBreakCellB = dict()
    rightGeneCellB = dict()
    rightGeneBreakCellB = dict()
    annotsCellB = dict()
    leftAndRightGenesB = list()

    #all necessary containers for comparison for cell line C if it exists
    fusionNameCellC = set()
    junctionReadCountCellC = dict()
    spanningFragCountCellC= dict()
    spliceTypeCellC = dict()
    leftGeneCellC = dict()
    leftGeneBreakCellC = dict()
    rightGeneCellC = dict()
    rightGeneBreakCellC = dict()
    annotsCellC = dict()
    leftAndRightGenesC = list()

    #all necessary containers for comparison for cell line D if it exists
    fusionNameCellD = set()
    junctionReadCountCellD = dict()
    spanningFragCountCellD= dict()
    spliceTypeCellD = dict()
    leftGeneCellD = dict()
    leftGeneBreakCellD = dict()
    rightGeneCellD = dict()
    rightGeneBreakCellD = dict()
    annotsCellD = dict()
    leftAndRightGenesD = list()


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
                junctionReadCountCellA.update({parts[0]: parts[1]})
                spanningFragsCellA.update({parts[0]: parts[2]})
                spliceTypeCellA.update({parts[0]: parts[3]})
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
                junctionReadCountCellB.update({parts[0]: parts[1]})
                spanningFragsCellB.update({parts[0]: parts[2]})
                spliceTypeCellB.update({parts[0]: parts[3]})
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
                    junctionReadCountCellC.update({parts[0]: parts[1]})
                    spanningFragsCellC.update({parts[0]: parts[2]})
                    spliceTypeCellC.update({parts[0]: parts[3]})
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
                    junctionReadCountCellD.update({parts[0]: parts[1]})
                    spanningFragsCellD.update({parts[0]: parts[2]})
                    spliceTypeCellD.update({parts[0]: parts[3]})
                    leftGeneCellD.update({parts[0]: left_rightGene[0]})
                    leftGeneBreakCellD.update({parts[0]: left_break[1]})
                    rightGeneCellD.update({parts[0]: left_rightGene[1]})
                    rightGeneBreakCellD.update({parts[0]: right_break[1]})
                    annotsCellD.update({parts[0]: parts[16]})
                    leftAndRightGenesD.extend((left_rightGene[0], left_rightGene[1]))

    #add to each set the fusions that are shared among the cell lines being compared
        AtoBNameConnections = intersection(fusionNameCellA, fusionNameCellB)
        AtoCNameConnections = intersection(fusionNameCellA, fusionNameCellC)
        AtoDNameConnections = intersection(fusionNameCellA, fusionNameCellD)
        BtoCNameConnections = intersection(fusionNameCellB, fusionNameCellC)
        BtoDNameConnections = intersection(fusionNameCellB, fusionNameCellD)
        CtoDNameConnections = intersection(fusionNameCellD, fusionNameCellD)
        AtoBtoCNameConnections = intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC)
        AtoBtoDNameConnections = intersection(fusionNameCellA, fusionNameCellB, fusionNameCellD)
        AtoCtoDNameConnections = intersection(fusionNameCellA, fusionNameCellC, fusionNameCellD)
        BtoCtoDNameConnections = intersection(fusionNameCellB, fusionNameCellC, fusionNameCellD)
        AtoBtoCtoDNameConnections = intersection(fusionNameCellA, fusionNameCellB, fusionNameCellC, fusionNameCellD)



    #create a Counter that will organize the genes into most frequent
    #the counter is a dictionary where each key is each gene and the value will be its frequency
        mostComGenesA = Counter(leftAndRightGenesA)
        mostComGenesB = Counter(leftAndRightGenesB)
        mostComGenesC = Counter(leftAndRightGenesC)
        mostComGenesD = Counter(leftAndRightGenesD)
        comGenesA_B = intersection(set(leftAndRightGenesA), set(leftAndRightGenesB))
        comGenesA_C = intersection(set(leftAndRightGenesA), set(leftAndRightGenesC))
        comGenesA_D = intersection(set(leftAndRightGenesA), set(leftAndRightGenesD))
        comGenesB_C = intersection(set(leftAndRightGenesB), set(leftAndRightGenesC))
        comGenesB_D = intersection(set(leftAndRightGenesB), set(leftAndRightGenesD))
        comGenesC_D = intersection(set(leftAndRightGenesC), set(leftAndRightGenesD))

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









if __name__ == '__main__':
    main()
