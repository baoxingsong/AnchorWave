#!python
import numpy as np
import sys
import re


#
# This script takes two known-splicesite-infile and novel-splicesite-outfile as input,
# try to figure out which known-splicesite-infile is more similar with novel-splicesite-outfile
# known-splicesite-infile was used as input for RNA-seq mapping using hisat2.
# novel-splicesite-outfile is the output of hisat2.
# Here I implemented this script to compare the output of function gffCoordinateLiftOver and the output of function annotationAndExonerateAndNovo
# For this aim I only consider the splice site of CDS elements
#


class Variant:
    def __init__(self, chr, position, ref, alt):
        self.chr=chr
        self.position=position
        self.ref=ref
        self.alt = alt
    def __hash__(self):
        return hash(self.chr) + hash(self.position) +  hash(self.ref) + hash(self.alt)

    def __lt__(self, other):
        if (self.chr < other.chr):
            return True
        if (self.chr == other.chr) and (self.position < other.position):
            return True
        if (self.chr == other.chr) and (self.position == other.position):
            return False
        return False

    def __gt__(self, other):
        if (self.chr > other.chr):
            return False
        if (self.chr == other.chr) and (self.position > other.position):
            return False
        if (self.chr == other.chr) and (self.position == other.position):
            return False
        return True

    def __eq__(self, other):
        return (self.chr == other.chr and self.position == other.position and self.ref == other.ref and self.alt == other.alt )


def readSdi(sdiFile):
    ss = set()
    with open(sdiFile) as f:
        for line in f:
            line = re.sub(r"\n", "", line)
            elements = line.split("\t")
            chr = elements[0]
            position = int(elements[1])
            ref=elements[3]
            alt = elements[4]

            ref = ref.replace( 'U', 'T')
            ref = ref.replace( 'R', 'N')
            ref = ref.replace( 'Y', 'N')
            ref = ref.replace( 'S', 'N')
            ref = ref.replace( 'W', 'N')
            ref = ref.replace( 'K', 'N')
            ref = ref.replace( 'M', 'N')
            ref = ref.replace( 'B', 'N')
            ref = ref.replace( 'D', 'N')
            ref = ref.replace( 'H', 'N')
            ref = ref.replace( 'V', 'N')

            alt = alt.replace( 'U', 'T')
            alt = alt.replace( 'R', 'N')
            alt = alt.replace( 'Y', 'N')
            alt = alt.replace( 'S', 'N')
            alt = alt.replace( 'W', 'N')
            alt = alt.replace( 'K', 'N')
            alt = alt.replace( 'M', 'N')
            alt = alt.replace( 'B', 'N')
            alt = alt.replace( 'D', 'N')
            alt = alt.replace( 'H', 'N')
            alt = alt.replace( 'V', 'N')
            ss.add( Variant( chr, position, ref, alt) )
    return ss

def readVcf(vcfFile):
    ss = set()
    with open(vcfFile) as f:
        for line in f:
            if line[0] != "#":
                line = re.sub(r"\n", "", line)
                elements = line.split("\t")
                chr = elements[0]
                position = int(elements[1])
                ref=elements[3]
                alt = elements[4]
                if len(ref) == 1 and len(alt)>1 and alt[0]==ref[0]:
                    position = position+1
                    ref = "-"
                    alt = alt[1:]
                elif len(ref) > 1 and len(alt) == 1 and alt[0]==ref[0]:
                    position = position+1
                    ref = ref[1:]
                    alt = "-"
                ss.add( Variant( chr, position, ref, alt) )
    return ss

ref = "ATWCRGRW"
ref = ref.replace( 'R', 'N')
ref = ref.replace( 'W', 'N')
print(ref)

sdi = readSdi("can_0.v7c.sdi")
totalLength = 0
for v in sdi:
    if v.ref != "-" :
        totalLength = totalLength + len(v.ref)

print(totalLength)


vcf = readVcf("can_0.vcf")
totalLength = 0
for v in vcf:
    if v.ref != "-" :
        totalLength = totalLength + len(v.ref)

print(totalLength)

totalLength = 0
for v in vcf:
    if (v not in sdi):
        if v.ref != "-" :
            totalLength = totalLength + len(v.ref)
print(totalLength)


for v in sdi:
    if (v not in vcf):
        print ( v.chr + "\t" + str(v.position) + "\t" + v.ref + "\t" + v.alt)



diffVariant = []
for v in vcf:
    if (v not in sdi) and (len(v.ref)>1 or len(v.alt)>1):
        v.alt = v.alt + "\tvcf"
        diffVariant.append(v)


for v in sdi:
    if (v not in vcf) and (len(v.ref)>1 or len(v.alt)>1):
        diffVariant.append(v)
        v.alt = v.alt + "\tsdi"


lastVcf = False
lastInsertion = False
diffVariant.sort()
for i in range(len(diffVariant)):
    v = diffVariant[i]
    if "vcf" in v.alt:
        if lastVcf:
            if lastInsertion and  len(v.ref) > len(v.alt):
                print ( diffVariant[i-1].chr + "\t" + str(diffVariant[i-1].position) + "\t" + diffVariant[i-1].ref + "\t" + diffVariant[i-1].alt)
                print ( v.chr + "\t" + str(v.position) + "\t" + v.ref + "\t" + v.alt)
            elif (not lastInsertion) and  len(v.ref) < len(v.alt):
                print ( diffVariant[i-1].chr + "\t" + str(diffVariant[i-1].position) + "\t" + diffVariant[i-1].ref + "\t" + diffVariant[i-1].alt)
                print ( v.chr + "\t" + str(v.position) + "\t" + v.ref + "\t" + v.alt)
        lastVcf = True
        if len(v.ref) > len(v.alt):
            lastInsertion = False
        else:
            lastInsertion = True
    else:
        lastVcf = False
