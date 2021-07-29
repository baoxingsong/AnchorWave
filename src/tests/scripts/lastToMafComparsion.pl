#!perl

use strict;
use warnings FATAL => 'all';
print "##maf version=1\n";
my $lineNUmber=0;
open(INPUT, '<',$ARGV[0]) or die $!;
while( my $line=<INPUT> ){
    if( $line=~/^#/ ){

    }else{
        if( $line=~/s\sChr/ ){
            if( $lineNUmber == 0 ){
                $line=~s/s\sChr/s\tcol.Chr/;
                $lineNUmber=1;
            }else{
                $line=~s/s\sChr/s\tquery.Chr/;
                $lineNUmber=0;
            }
        }
        print "$line";
    }
}
close INPUT;
