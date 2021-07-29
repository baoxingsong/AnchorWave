#!perl -w
use strict;
open INPUT, "$ARGV[1]";
my $totalNUmber = 0;
my $goodnumber = 0;
while( my $line=<INPUT> ){
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	my @words = split /\s+/, $line;
	my $chr = $words[0];
	my $start = $words[1];
	$start = $start + 1;
	my $end = $words[2];
	my $teLength = $end - $start + 1;

	$start = $start - 35;
	$end = $end + 35;

	my $maximumGapLength = 0;
	my $maximumGapPreNonGapLength = 0;
	my $maximumGapPostNonGapLength = -1;

	my $thisGapLength = 0;
	my $thisNonGapLength = 0;
	my $lastGap = 0;

	#my $result = `samtools depth -r $chr:$start-$end $ARGV[0]`;
	my $result = `samtools mpileup -r $chr:$start-$end $ARGV[0]`;
	my @lines = split /\n+/, $result;

	my $lastPosition = $start - 1;
	
	my $numberOfLines = scalar(@lines);
	my $nonZoroCoverage = 0;
	for my $eachline (@lines){
		my @elements = split /\s+/, $eachline;
		# for( my $i= $lastPosition+1; $i<$elements[1]; $i = $i + 1){
		# 	if( $lastGap ){
		# 		$thisGapLength = $thisGapLength + 1;
		# 	}else{
		# 		if( $maximumGapPostNonGapLength == -1 ){
		# 			$maximumGapPostNonGapLength = $thisNonGapLength;
		# 		}
		# 		$thisGapLength = 1;
		# 	}
		# 	$lastGap = 1;
		# 	$numberOfLines = $numberOfLines + 1;
		# }
		$lastPosition = $elements[1];

		#if ( $elements[2] != 0 ){
		if ( $elements[4] ne "*" ){
			$nonZoroCoverage = $nonZoroCoverage + 1;
			
			if( $lastGap ){
				if( $thisGapLength > $maximumGapLength ){
					$maximumGapLength = $thisGapLength;
					$maximumGapPreNonGapLength = $thisNonGapLength;
					$maximumGapPostNonGapLength = -1;
				}
				$thisNonGapLength = 1;
			}else{
				$thisNonGapLength = $thisNonGapLength + 1;
				
			}
			$lastGap = 0;
		}else{
			if( $lastGap ){
				$thisGapLength = $thisGapLength + 1;
			}else{
				if( $maximumGapPostNonGapLength == -1 ){
					$maximumGapPostNonGapLength = $thisNonGapLength;
				}
				$thisGapLength = 1;
			}
			$lastGap = 1;
		}
	}
	if( $maximumGapPostNonGapLength == -1 ){
		$maximumGapPostNonGapLength = $thisNonGapLength;
	}
	if( $thisGapLength > $maximumGapLength ){
		$maximumGapLength = $thisGapLength;
		$maximumGapPostNonGapLength = 0;
	}

	print "$line\t$teLength\t$maximumGapLength\t$maximumGapPreNonGapLength\t$maximumGapPostNonGapLength\t";

	if( ($teLength + 70) > $numberOfLines ){
		print "IncompleteAlignment\n";
		$totalNUmber = $totalNUmber + 1;
	}elsif( $maximumGapPostNonGapLength == 0 || $maximumGapPreNonGapLength == 0 ){
		print "BOUNDARIESNOTCovered\n";
		$totalNUmber = $totalNUmber + 1;
	}elsif ( $maximumGapLength < $teLength && ($teLength - $maximumGapLength) < 35 ){
		print "GOOD1\n";
		$totalNUmber = $totalNUmber + 1;
		$goodnumber = $goodnumber + 1;
	}elsif ( $maximumGapLength < $teLength ){
		print "TEBREAKED\n";
		$totalNUmber = $totalNUmber + 1;
	}elsif ( $maximumGapLength > $teLength && ($maximumGapLength - $teLength) < 35 ){
		print "GOOD2\n";
		$totalNUmber = $totalNUmber + 1;
		$goodnumber = $goodnumber + 1;
	}elsif ( $maximumGapLength > $teLength ){
		print "TooLargeIndeL\n";
		$totalNUmber = $totalNUmber + 1;
	}elsif ( $maximumGapPostNonGapLength > 0 && $maximumGapPreNonGapLength > 0 ){
		print "GOOD\n";
		$totalNUmber = $totalNUmber + 1;
		$goodnumber = $goodnumber + 1;
	}else{
		print "TOBECHECKED\n";
		$totalNUmber = $totalNUmber + 1;
	}
}
close INPUT;

print "totalNUmber:$totalNUmber\n";
print "goodnumber:$goodnumber\n";
