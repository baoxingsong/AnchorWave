#!python3
import sys
import re
import os

class samRecord:
    def __init__(self, wholeLine, chr, start, end):
        self.wholeLine = wholeLine
        self.chr = chr
        self.start = start
        self.end = end


def overLap(samRecord1, samRecord2):
    if samRecord1.chr == samRecord2.chr:
        if samRecord2.start <= samRecord1.start and samRecord1.start <= samRecord2.end:
            return True
        if samRecord2.start <= samRecord1.end and samRecord1.end <= samRecord2.end:
            return True
        if samRecord1.start <= samRecord2.start and samRecord2.start <= samRecord1.end:
            return True
        if samRecord1.start <= samRecord2.end and samRecord2.end <= samRecord1.end:
            return True
    return False

def readsam(sam_file):
    samRecords = []
    with open(sam_file) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search('^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s', line)
                if( m != None and m.group(3) != "*" ):
                    chromosome_name = m.group(3)
                    start = int(m.group(4))
                    end = int(m.group(4))
                    ms = re.findall('(\d+)[MD=X]', m.group(6))
                    for mm in ms:
                        end += int(mm)
                    if( end > start ):
                        end -= 1
                    samRecords.append(samRecord(line, chromosome_name, start, end))
    return samRecords


samRecords = readsam(sys.argv[1])
print (str(len(samRecords)))
id = 0
while len(samRecords) > 0:
    toOutPut = []
    toKept = []
    i = 0
    toOutPut.append(samRecords[i])
    for i in range(1, len(samRecords)):
        overlapped = False
        for samRecord1 in toOutPut:
            if overLap(samRecord1, samRecords[i]):
                overlapped = True
                break
        if overlapped:
            toKept.append(samRecords[i])
        else:
            toOutPut.append(samRecords[i])
    samFile = sys.argv[1] + "." + str(id) + ".sam"
    bamFile = sys.argv[1] + "." + str(id) + ".bam"
    output = open(samFile, 'w')
    for samRecord1 in toOutPut:
        output.write(samRecord1.wholeLine)
    output.close()
    command = "samtools view -O BAM --reference Zea_mays.AGPv4.dna.toplevel.fa " + samFile + " | samtools sort - > " + bamFile
    os.system(command)
    id = id + 1
    samRecords = toKept
