#!perl
use strict;
use warnings FATAL => 'all';

my $vcfFile = $ARGV[0];

open INPUT, "$vcfFile";
while( my $line=<INPUT> ){
    $line=~s/^chr//g;
    if( substr($line, 0, 1) ne "#" ){
        my @element = split('\t', $line);
        my $length = length($element[4]) - length($element[3]);
        my $chr = $element[0];
        my $position = $element[1];
        my $ref = $element[3];
        my $alt = $element[4];
#        print "$chr\t$position\t$length\t$ref\t$alt\n";
        if( substr($ref, 0, 1) eq substr($alt, 0, 1) ){
            $position = $position - 1;
            $ref = substr $ref, 1;
            $alt = substr $alt, 1;
            if(length($ref)==0){
                $ref = "-";
            }
            if(length($alt)==0){
                $alt = "-";
            }
 #           print "$chr\t$position\t$length\t$ref\t$alt\n";
        }elsif( length($ref)==1 && length($alt)==1 ){
        }elsif( (substr $ref, -1) eq (substr $alt, -1) ){
            #print "$line";
            chop($ref);
            chop($alt);
            if(length($ref)==0){
                $ref = "-";
            }
            if(length($alt)==0){
                $alt = "-";
            }
            # print "$chr\t$position\t$length\t$ref\t$alt\n";
        }else{
#            print "$line";
        }
        print "$chr\t$position\t$length\t$ref\t$alt\n";
    }
}
close INPUT;

