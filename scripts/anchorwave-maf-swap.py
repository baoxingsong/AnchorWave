#! /usr/bin/env python

# Read MAF-format alignments, and write them, after moving the Nth
# sequence to the top in each alignment.

# Before writing, if the top sequence would be on the - strand, then
# flip all the strands.  But don't do this if the top sequence is
# translated DNA.

from __future__ import print_function

import gzip
import itertools
import optparse
import os
import signal
import sys

try:
    maketrans = str.maketrans
except AttributeError:
    from string import maketrans

def myOpen(fileName):
    if fileName == "-":
        return sys.stdin
    if fileName.endswith(".gz"):
        return gzip.open(fileName, "rt")  # xxx dubious for Python2
    return open(fileName)

def indexOfNthSequence(mafLines, n):
    for i, line in enumerate(mafLines):
        if line[0] == "s":
            if n == 1: return i
            n -= 1
    raise Exception("encountered an alignment with too few sequences")

def rangeOfNthSequence(mafLines, n):
    """Get the range of lines associated with the Nth sequence."""
    start = indexOfNthSequence(mafLines, n)
    stop = start + 1
    while stop < len(mafLines):
        line = mafLines[stop]
        if line[0] not in "qi": break
        stop += 1
    return start, stop

complement = maketrans('ACGTNSWRYKMBDHVacgtnswrykmbdhv',
                       'TGCANSWYRMKVHDBtgcanswyrmkvhdb')
# doesn't handle "U" in RNA sequences
def revcomp(seq):
    return seq[::-1].translate(complement)

def flippedMafRecords(mafLines):
    for line in mafLines:
        words = line.split()
        if words[0] == "s":
            s, name, start, span, strand, seqlen, alnString = words
            #newStart = str(int(seqlen) - int(start) - int(span))
            newStart = start
            newStrand = "+-"[strand == "+"]
            newString = revcomp(alnString)
            yield [s, name, newStart, span, newStrand, seqlen, newString]
        elif words[0] == "p":
            yield words[:1] + [words[1][::-1]]
        elif words[0] == "q":
            yield words[:2] + [words[2][::-1]]
        else:
            yield words

def sLineFieldWidths(mafLines):
    sLines = (i for i in mafLines if i[0] == "s")
    sColumns = zip(*sLines)
    for i in sColumns:
        yield max(map(len, i))

def joinedMafS(fieldWidths, words):
    formatParams = itertools.chain.from_iterable(zip(fieldWidths, words))
    return "%*s %-*s %*s %*s %*s %*s %*s\n" % tuple(formatParams)

def joinedMafLine(words, fieldWidths):
    if words[0] == "s":
        return joinedMafS(fieldWidths, words)
    elif words[0] == "q":
        words = words[:2] + [""] * 4 + words[2:]
        return joinedMafS(fieldWidths, words)
    elif words[0] == "p":
        words = words[:1] + [""] * 5 + words[1:]
        return joinedMafS(fieldWidths, words)
    else:
        return " ".join(words) + "\n"

def flippedMaf(mafLines):
    flippedLines = list(flippedMafRecords(mafLines))
    fieldWidths = list(sLineFieldWidths(flippedLines))
    return (joinedMafLine(i, fieldWidths) for i in flippedLines)

def isCanonicalStrand(mafLine):
    words = mafLine.split()
    strand = words[4]
    if strand == "+": return True
    alnString = words[6]
    if "/" in alnString or "\\" in alnString: return True  # frameshifts
    alnSize = int(words[3])
    gapCount = alnString.count("-")
    if len(alnString) - gapCount < alnSize: return True  # translated DNA
    return False

def swapOneMaf(opts, mafLines):
    start, stop = rangeOfNthSequence(mafLines, opts.n)
    mafLines[1:stop] = mafLines[start:stop] + mafLines[1:start]
    if not isCanonicalStrand(mafLines[1]):
        mafLines = flippedMaf(mafLines)
    for i in mafLines:
        print(i, end="")
    print()  # blank line after each alignment

def mafSwap(opts, args):
    if not args:
        args = ["-"]
    for fileName in args:
        inputLines = myOpen(fileName)
        mafLines = []
        for line in inputLines:
            if line[0] == "#":
                print(line, end="")
            elif line.isspace():
                if mafLines:
                    swapOneMaf(opts, mafLines)
                    mafLines = []
            else:
                mafLines.append(line)
        if mafLines:
            swapOneMaf(opts, mafLines)

if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # avoid silly error message

    usage = "%prog [options] my-alignments.maf"
    description = "Change the order of sequences in MAF-format alignments."
    op = optparse.OptionParser(usage=usage, description=description)
    op.add_option("-n", type="int", default=2,
                  help="move the Nth sequence to the top (default: %default)")
    (opts, args) = op.parse_args()
    if opts.n < 1: op.error("option -n: should be >= 1")

    try: mafSwap(opts, args)
    except KeyboardInterrupt: pass  # avoid silly error message
    except Exception as e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
