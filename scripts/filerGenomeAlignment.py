#!python
import numpy as np
import sys
import re

class Cigar:
    def __init__(self, length, type, refPosition, quePosition):
        self.length = length
        self.type = type
        self.refPosition = refPosition
        self.quePosition = quePosition


def findGood(cigars, queryName, flag, referenceName, mapq, queryPosition, cigar, queryPositionEnd, queSeq, matchScore, mismatchScore, scoreThreshold, xscore):
    maximumScore = 0
    thisScore = 0
    startIndex = 0
    endIndex = 0
    currentQuerySeqPosition = 0
    endQuerySeqPosition = 0
    startQuerySeqPosition = 0

    i = 0
    while i < len(cigars):
        cigar = cigars[i]
        if cigar.type == "=":
            thisScore = thisScore + matchScore * cigar.length
            currentQuerySeqPosition = currentQuerySeqPosition + cigar.length
        elif cigar.type == "X":
            thisScore = thisScore - mismatchScore * cigar.length
            currentQuerySeqPosition = currentQuerySeqPosition + cigar.length
        else:
            thisScore = thisScore - mismatchScore
            if cigar.type == "I":
                currentQuerySeqPosition = currentQuerySeqPosition + cigar.length
        if maximumScore < thisScore:
            maximumScore = thisScore
            endIndex = i
            endQuerySeqPosition = currentQuerySeqPosition
        i = i + 1
        # if (maximumScore-thisScore) >= xscore:
        #     if maximumScore >= scoreThreshold:
        #         i = len(cigars)
        #     else:
        #         thisScore = 0
        #         maximumScore = 0
        if thisScore < 0:
            thisScore = 0

    currentQuerySeqPosition = endQuerySeqPosition
    # print("maximumScore:" + str(maximumScore))
    # print ("startIndex:" + str(startIndex))
    # print ("endIndex:" + str(endIndex))
    # print (queryName + "\t" + str(flag) + "\t" + referenceName + "\t" + str(cigars[startIndex].refPosition) + "\t" + str(mapq) + "\t", end = '')
    if maximumScore >= scoreThreshold:
        maximumScore = 0
        thisScore = 0
        i = endIndex
        while i >= 0:
            cigar = cigars[i]
            if cigar.type == "=":
                thisScore = thisScore + matchScore * cigar.length
                currentQuerySeqPosition = currentQuerySeqPosition - cigar.length
            elif cigar.type == "X":
                thisScore = thisScore - mismatchScore * cigar.length
                currentQuerySeqPosition = currentQuerySeqPosition - cigar.length
            else:
                thisScore = thisScore - mismatchScore
                if cigar.type == "I":
                    currentQuerySeqPosition = currentQuerySeqPosition - cigar.length

            if maximumScore < thisScore:
                maximumScore = thisScore
                startIndex = i
                startQuerySeqPosition = currentQuerySeqPosition
            i = i - 1
            # if (maximumScore-thisScore) >= xscore:
            #     i = -1
        if maximumScore >= scoreThreshold:
            # print ("startIndex:" + str(startIndex))
            # print("maximumScore:" + str(maximumScore))
            print (queryName + "\t" + str(flag) + "\t" + referenceName + "\t" + str(cigars[startIndex].refPosition) + "\t" + str(mapq) + "\t", end = '')
            i = startIndex
            while i <= endIndex:
                print (str(cigars[i].length) + cigars[i].type, end = '')
                i = i + 1
            print ("\t*\t0\t0\t" + queSeq[startQuerySeqPosition:endQuerySeqPosition] + "\t*")
            #print ("\t*\t0\t0\t*" + "\t*")
            if startIndex > 0:
                findGood(cigars[0:startIndex], queryName, flag, referenceName, mapq, queryPosition, cigar, queryPositionEnd, queSeq[0:startQuerySeqPosition], matchScore, mismatchScore, scoreThreshold, xscore)
            if endIndex+1 < len(cigars):
                findGood(cigars[endIndex+1:len(cigars)], queryName, flag, referenceName, mapq, queryPosition, cigar, queryPositionEnd, queSeq[endQuerySeqPosition:len(queSeq)], matchScore, mismatchScore, scoreThreshold, xscore)
    return 0

def parseSam(queryName, flag, referenceName, referencePosition, mapq, queryPosition, cigar, queryPositionEnd, queSeq, matchScore=37, mismatchScore=63, scoreThreshold=999, xscore=54):
    pattern = re.compile(r'([0-9]+)([A-Z=]+)')
    thisRefPosition = referencePosition
    thisQuePosition = queryPosition
    if flag == 16:
        thisQuePosition = queryPositionEnd

    cigars = []
    for (numbers, type) in re.findall(pattern, cigar):
        numbers = int(numbers)
        cigars.append( Cigar( numbers, type, thisRefPosition, thisQuePosition))
        if flag == 0 :
            if type == "=" or type == "X" :
                thisRefPosition = thisRefPosition + numbers
                thisQuePosition = thisQuePosition + numbers
            elif type == "D":
                thisRefPosition = thisRefPosition + numbers
            elif type == "I":
                thisQuePosition = thisQuePosition + numbers
            else:
                print ("unsupported cigar type:" + type)
        elif flag == 16:
            if type == "=" or type == "X" :
                thisRefPosition = thisRefPosition + numbers
                thisQuePosition = thisQuePosition - numbers
            elif type == "D":
                thisRefPosition = thisRefPosition + numbers
            elif type == "I":
                thisQuePosition = thisQuePosition - numbers
            else:
                print ("unsupported cigar type:" + type)
        else:
            print ("unsupported flag:" + flag)

    if len(cigars) == 1:
        for cig in cigars:
            print (queryName + "\t" + str(flag) + "\t" + referenceName + "\t" + str(referencePosition) + "\t" + str(mapq) + "\t" + str(queryPosition) + "H" + str(cig.length) + cig.type + "\t*\t0\t0\t" + queSeq + "\t*")
    else:
        # cig = cigars[0]
        # print(cig.refPosition)
        # cig = cigars[1]
        # print(cig.refPosition)
        findGood(cigars, queryName, flag, referenceName, mapq, queryPosition, cigar, queryPositionEnd, queSeq, matchScore, mismatchScore, scoreThreshold, xscore)

with open("/media/bs674/ppi8t/testWAF/alignSorghumAgainstMaizeV4/newversions/lm23v2/alignment.sam") as f:
    for line in f:
        elements = re.match(r'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)H(\S+)(\d+)H\s+\S+\s+\S+\s+\S+\s+(\S+)\s', line)
        if elements:
            parseSam(elements.group(1), int(elements.group(2)), elements.group(3), int(elements.group(4)), int(elements.group(5)), int(elements.group(6)), elements.group(7), int(elements.group(8)), elements.group(9))
