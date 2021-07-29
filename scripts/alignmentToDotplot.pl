#!perl
use strict;
use warnings FATAL => 'all';

my $transcript_gff_file = $ARGV[0];
my $sam_file = $ARGV[1];

my %transcriptChr;
my %transcriptPosition;
my %transcriptStrand;
open INPUT, "$transcript_gff_file";
while( my $line=<INPUT> ){
	if( $line=~/^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+\d+\s+\S+\s+(\S+)\s+\S+\s.*ID=([\.a-zA-Z0-9\:\-_]*?)(;|\z)/ ){
		$transcriptChr{$5}=$1;
		$transcriptPosition{$5}=$3;
		$transcriptStrand{$5}=$4;
	}
}
close INPUT;

open INPUT, "$sam_file";
while( my $line=<INPUT> ){
	if( $line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+/ ){
	    if( exists $transcriptStrand{$1} ){
            if( ( $2%32 ==0 && ($transcriptStrand{$1} eq "+") ) || ( $2%32 ==16 && ($transcriptStrand{$1} eq "-") ) ){
                print "$transcriptChr{$1}\t$transcriptPosition{$1}\t$3\t$4\t+\n";
            }else{
                print "$transcriptChr{$1}\t$transcriptPosition{$1}\t$3\t$4\t-\n";
            }
         }else{
            print STDERR "could not find $1\n";
         }
	}
}
close INPUT;
