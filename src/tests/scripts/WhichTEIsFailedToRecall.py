#!python3
import sys

class Anchor:
    def __init__(self, refChr, refStart, refEnd, queryChr, queryStart, queryEnd, name):
        self.refChr = refChr
        self.refStart = refStart
        self.refEnd = refEnd
        self.queryChr = queryChr
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.name = name

class Deletions:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end
        self.overlapped = ""

def readAnchors(anchorFile):
    anchors = dict()
    with open(anchorFile) as f:
        linenumber = 0
        for line in f:
            if line[0] != "#" and line[:3] != "ref":
                (refChr, refStart, refEnd, queryChr, queryStart, queryEnd, strand, name, blockid, similarity) = line.split()
                if refChr not in anchors:
                    anchors[refChr] = []
                refStart = int(refStart)
                refEnd = int(refEnd)
                if 0 == linenumber:
                    queryStart = int(queryStart)
                    anchors[refChr].append(Anchor(refChr, 1, refStart-1, queryChr, 1, queryStart-1, name))
                    print("line 35")
                anchors[refChr].append(Anchor(refChr, refStart, refEnd, queryChr, queryStart, queryEnd, name))
                linenumber = linenumber + 1
    return anchors

def readErrors(logFile):
    errorDeletioins = []
    with open(logFile) as f:
        for line in f:
            (code, chr, position, length) = line.split()[:4]
            position = int(position)
            length = -int(length)
            if code == "bad":
                errorDeletioins.append(Deletions(chr, position, position+length-1))
            elif code != "good":
                print("hahah" + line)
    return errorDeletioins



def readCorrect(logFile):
    errorDeletioins = []
    with open(logFile) as f:
        for line in f:
            (code, chr, position, length) = line.split()[:4]
            position = int(position)
            length = -int(length)
            if code == "good":
                errorDeletioins.append(Deletions(chr, position, position+length-1))
    return errorDeletioins

fa = 80000
fa2 = 120000
anchors = readAnchors(sys.argv[1])
errorDeletioins = readErrors(sys.argv[2])
correctDeletioins = readCorrect(sys.argv[2])
print (len(errorDeletioins))
print (len(correctDeletioins))
for deletion in errorDeletioins:
    numberOfOverlappedAnchors = 0
    theProblemAnchor = "NA"
    for anchor in anchors[deletion.chr]:
        if (anchor.refStart <= deletion.start and deletion.start <= anchor.refEnd) or (anchor.refStart <= deletion.end and deletion.end <= anchor.refEnd) or (deletion.start<=anchor.refStart and anchor.refStart<=deletion.end):
            deletion.overlapped = deletion.overlapped + "_" + anchor.name
            numberOfOverlappedAnchors = numberOfOverlappedAnchors + 1
            theProblemAnchor = anchor
            #print (deletion.chr + "\t" + str(deletion.start) + "\t" + str(deletion.end) + "\t" + anchor.refChr + "\t" + str(anchor.refStart) + "\t" + str(anchor.refEnd) + "\t" + anchor.queryChr + "\t" + anchor.queryStart + "\t" + anchor.queryEnd + "\t" + anchor.name)
    if numberOfOverlappedAnchors > 1:
        print ("bad\t" + deletion.chr + "\t" + str(deletion.start) + "\t" + str(deletion.end) + "\tboundary\t" +  deletion.overlapped +"\t0")
    else:
        anchorLength = theProblemAnchor.refEnd -  theProblemAnchor.refStart + 1
        if "interanchor" == theProblemAnchor.name and anchorLength > fa:
            method = "slidingWindow"
        elif "interanchor" != theProblemAnchor.name and anchorLength > fa2:
            method = "slidingWindow"
        else:
            method = "WFA"
        print ("bad\t" + deletion.chr + "\t" + str(deletion.start) + "\t" + str(deletion.end) + "\t" + method + "\t" + theProblemAnchor.name+"\t" + str(theProblemAnchor.refEnd- theProblemAnchor.refStart+1))

for deletion in correctDeletioins:
    overlappedNUmber = 0
    for anchor in anchors[deletion.chr]:
        if (anchor.refStart <= deletion.start and deletion.start <= anchor.refEnd) or (anchor.refStart <= deletion.end and deletion.end <= anchor.refEnd):
            overlappedNUmber = overlappedNUmber + 1
            anchorLength = anchor.refEnd -  anchor.refStart + 1
            if "interanchor" == anchor.name and anchorLength > fa:
                method = "slidingWindow"
            elif "interanchor" != anchor.name and anchorLength > fa2:
                method = "slidingWindow"
            else:
                method = "WFA"
            print ("good\t" + deletion.chr + "\t" + str(deletion.start) + "\t" + str(deletion.end) + "\t" + method + "\t" + anchor.name+"\t" + str(anchor.refEnd- anchor.refStart+1))
