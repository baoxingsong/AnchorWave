#!python3
import sys

def readMaf(benchMark_file, test_file):
    benchMarkalignedPosition = dict()
    with open(benchMark_file) as f:
        for line in f:
            elements = line.split()
            if len(elements) == 7:
                (s, refchr, refstart, reflength, refstrand, refchrLength, refali) = elements
                refchr = refchr.replace("col.", "")
                if refchr not in benchMarkalignedPosition:
                    benchMarkalignedPosition[refchr] = dict()
                refstart = int(refstart)

                line2 = f.readline()
                elements2 = line2.split()
                (s, querychr, querystart, querylength, querystrand, queryChrLength, queryali) = elements2
                querychr = querychr.replace("query.", "")
                if querychr != refchr:
                    print("this script does not apply for you comparsion, please do not use it")
                    print(querychr)
                    print(refchr)
                querystart = int(querystart)
                refPosition = 0
                queryPosition = 0
                if querystrand[0] == '+':
                    for i in range(len(refali)):
                        if refali[i] != '-':
                            if queryali[i] != '-' and refali[i] != queryali[i]:
                                benchMarkalignedPosition[refchr][refstart+refPosition] = querystart + queryPosition
                            elif queryali[i] == '-':
                                benchMarkalignedPosition[refchr][refstart+refPosition] = -1 # gap alignment, give it a value of -1
                        if refali[i] != '-':
                            refPosition = refPosition + 1
                        if queryali[i] != '-':
                            queryPosition = queryPosition + 1
                else:
                    print("this script does not apply for you comparsion, please do not use it")

    totalProducedLength = 0
    with open(test_file) as f:
        for line in f:
            elements = line.split()
            if len(elements) == 7:
                (s, refchr, refstart, reflength, refstrand, refchrLength, refali) = elements
                refchr = refchr.replace("col.", "")
                if refchr in benchMarkalignedPosition:
                    refstart = int(refstart)

                    line2 = f.readline()
                    elements2 = line2.split()
                    (s, querychr, querystart, querylength, querystrand, queryChrLength, queryali) = elements2
                    if len(refali) == len(queryali):
                        querychr = querychr.replace("query.", "")
                        if querychr != refchr:   ## interchrosome relocations are all wrong
                            for i in range(len(refali)):
                                if refali[i] != '-':
                                    totalProducedLength = totalProducedLength + 1
                        else:
                            querystart = int(querystart)
                            refPosition = 0
                            queryPosition = 0
                            if querystrand[0] == '+':
                                for i in range(len(refali)):
                                    if refali[i] != '-':
                                        if queryali[i] != '-' and (refstart+refPosition) in benchMarkalignedPosition[refchr] and benchMarkalignedPosition[refchr][refstart+refPosition] == (querystart + queryPosition):
                                            benchMarkalignedPosition[refchr][refstart+refPosition] = -5 # if it is same with benchmark, give it a value of -5
                                        elif queryali[i] == '-' and (refstart+refPosition) in benchMarkalignedPosition[refchr] and benchMarkalignedPosition[refchr][refstart+refPosition] == -1:
                                            benchMarkalignedPosition[refchr][refstart+refPosition] = -5
                                        if refali[i] != queryali[i]:
                                            totalProducedLength = totalProducedLength + 1
                                    if refali[i] != '-':
                                        refPosition = refPosition + 1
                                    if queryali[i] != '-':
                                        queryPosition = queryPosition + 1
                            else:    # here we code assume that there is no negative alignment in the benchmark alignment. All inversions are wrong If this is not true, should changed it
                                for i in range(len(refali)):
                                    if refali[i] != '-':
                                        totalProducedLength = totalProducedLength + 1
    return benchMarkalignedPosition, totalProducedLength


benchMarkalignedPosition, totalProducedLength = readMaf(sys.argv[1], sys.argv[2])

totalRefLength = 0
totalCorrected = 0
for chr in benchMarkalignedPosition:
    for position in benchMarkalignedPosition[chr]:
        totalRefLength = totalRefLength + 1
        if benchMarkalignedPosition[chr][position] == -5: # if it was same with benchmark, the value was set as -5
            totalCorrected = totalCorrected + 1

output = open(sys.argv[2] + ".aliEvaluatioin", 'w')
output.write ("totalRefLength:" + str(totalRefLength) + "\n")
output.write ("totalProducedLength:" + str(totalProducedLength) + "\n")
output.write ("totalCorrected:" + str(totalCorrected) + "\n")
recall = totalCorrected/totalRefLength
output.write ("recall:" + str(recall) + "\n")
precision = totalCorrected/totalProducedLength
output.write ("precision:" + str(precision) + "\n")
fscore = 2 * precision * recall/(precision + recall)
output.write ("fscore:" + str(fscore) + "\n")
output.close()
