#!perl

use strict;
use warnings FATAL => 'all';

print "# LAST version 1110\n";
print "#\n";
print "# a=7 b=1 A=7 B=1 e=30 d=21 x=29 y=9 z=29 D=1e+06 E=1760.84\n";
print "# R=10 u=0 s=2 S=0 M=0 T=0 m=10 l=1 n=10 k=1 w=1000 t=0.910239 j=3 Q=0\n";
print "# col\n";
print "# Reference sequences=5 normal letters=118960141\n";
print "# lambda=1.09602 K=0.335388\n";
print "#\n";
print "#     A   C   G   T   M   S   K   W   R   Y   B   D   H   V\n";
print "# A   1  -1  -1  -1   0  -1  -1   0   0  -1  -1   0   0   0\n";
print "# C  -1   1  -1  -1   0   0  -1  -1  -1   0   0  -1   0   0\n";
print "# G  -1  -1   1  -1  -1   0   0  -1   0  -1   0   0  -1   0\n";
print "# T  -1  -1  -1   1  -1  -1   0   0  -1   0   0   0   0  -1\n";
print "# M   0   0  -1  -1   0   0  -1   0   0   0   0   0   0   0\n";
print "# S  -1   0   0  -1   0   0   0  -1   0   0   0   0   0   0\n";
print "# K  -1  -1   0   0  -1   0   0   0   0   0   0   0   0   0\n";
print "# W   0  -1  -1   0   0  -1   0   0   0   0   0   0   0   0\n";
print "# R   0  -1   0  -1   0   0   0   0   0  -1   0   0   0   0\n";
print "# Y  -1   0  -1   0   0   0   0   0  -1   0   0   0   0   0\n";
print "# B  -1   0   0   0   0   0   0   0   0   0   0   0   0   0\n";
print "# D   0  -1   0   0   0   0   0   0   0   0   0   0   0   0\n";
print "# H   0   0  -1   0   0   0   0   0   0   0   0   0   0   0\n";
print "# V   0   0   0  -1   0   0   0   0   0   0   0   0   0   0\n";
print "#\n";
print "# Coordinates are 0-based.  For - strand matches, coordinates\n";
print "# in the reverse complement of the 2nd sequence are used.\n";
print "#\n";
print "# name start alnSize strand seqSize alignment\n";
print "#\n";
print "# batch 0\n";

my $lineNUmber=0;
open(INPUT, '<',$ARGV[0]) or die $!;
while( my $line=<INPUT> ){
    if( $line=~/^#/ ){

    }else{
        print "$line";
    }
}
close INPUT;
