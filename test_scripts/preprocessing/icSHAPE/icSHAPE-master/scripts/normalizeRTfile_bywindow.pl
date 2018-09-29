#! /usr/bin/env perl
#
#
use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use vars qw ($opt_h $opt_V $opt_D $opt_i $opt_o $opt_m $opt_d $opt_l $opt_f $opt_w $opt_s );
&getopts('hVDi:o:r:m:d:l:f:s:w:');

my $usage = <<_EOH_;
## --------------------------------------
Normalize base density file

Command:
$0 -i input_baseDensity_file -o output_normalized_baseDensity_file

# what it is:
 -i     input_baseDensity_file
 -o     output_normalized_baseDensity_file

# more options:
 -d     head to skip
 -l     tail to skip
 -f     scaling form (the benchmarked value will be scaled to this)
 -w     normalized windowsize (default: 200nt)
 -m     normalize method
 -s     step length, set it 0 to use disjoint window not sliding window (default: 50nt sliding window) 
_EOH_
;

&main();

sub main {
    my %parameters = &init();
    &normalizeBaseDensity ( $parameters{input}, $parameters{output}, $parameters{normalizeMethod}, $parameters{headToSkip}, $parameters{tailToSkip}, $parameters{scalingForm}, 
        $parameters{windowSize}, $parameters{stepSize} );
}

sub init {
    my %parameters = ();

    die $usage if ( $opt_h || ( not defined $opt_i ) || ( not defined $opt_o ) );
    $opt_V = 0 if ( not defined $opt_V );
    $opt_D = 0 if ( not defined $opt_D );

    $parameters{input} = $opt_i;
    $parameters{output} = $opt_o;

    if ( defined $opt_m ) { $parameters{normalizeMethod} = $opt_m; }
    else { $parameters{normalizeMethod} = "median:decile2"; }
    if ( defined $opt_d ) { $parameters{headToSkip} = $opt_d; }
    else { $parameters{headToSkip} = 0; }
    if ( defined $opt_l ) { $parameters{tailToSkip} = $opt_l; }
    else { $parameters{tailToSkip} = 0; }
    if ( defined $opt_f ) { $parameters{scalingForm} = $opt_f; }
    else { $parameters{scalingForm} = 100; }
    if ( defined $opt_w ) { $parameters{windowSize} = $opt_w; }
    else { $parameters{windowSize} = 200 ;}
    if ( defined $opt_s ) { 
        if ($opt_s == 0){
            $parameters{stepSize} = $parameters{windowSize} ;
        }
        else{
            $parameters{stepSize} = $opt_s; 
        }
        
    } 
    else { $parameters{stepSize} = 50; }

    return ( %parameters );
}


sub normalizeBaseDensity {
    my ( $baseDensityFile, $outputFile, $normalizeMethod, $headToSkip, $tailToSkip, $scalingForm, $windowSize, $stepSize ) = @_;
    my $outputFile_window = $outputFile.".win$windowSize.stp$stepSize" ;

    print STDERR "Normalize base density from file $baseDensityFile...\n\t", `date`;

    $headToSkip++; 
    open ( IN, $baseDensityFile );
    #open ( OUT, ">$outputFile" ) or die "Error opening $outputFile for writing transcript base density.\n";
    open ( OUT2, ">$outputFile_window" ) or die "Error opening $outputFile_window for writing transcript base density.\n";
    #print OUT "#transcript\tlength\ttype\tbase_frequency, start from position 1.\n";
    print OUT2 "#transcript\tlength\ttype\tbase_frequency, start from position 1. normalized by $windowSize window\n";
    
    #my %finalBaseDensity = ();
    #my %finalRTstop = ();
    my $transcript = "";  my $len = 0;  my $rpkm = 0; my @baseDensity = ();
    my $transCount = 0;
    while ( my $line = <IN> )  {
        next if ( $line =~ /^#/ );
        $transCount++ ;
        if (( $transCount % 5000) == 0 ){
            print STDERR "\tTranscript processed: $transCount\n\t", `date`;
        }
        my $scalingFactor = 1;
        chomp $line;
        ## remember that the first element in @baseDensity is base zero
        #  in RT stop file, it means no modified base in the read
        #  in background file, it means null, so should always be zero
        ( $transcript, $len, $rpkm, @baseDensity ) = split ( /\t/, $line );
        my $trimed_last = $len - $tailToSkip;      ## again, to avoid base zero
        while ( $trimed_last < ( $headToSkip + 40 ) ) {
            print STDERR "Warning! Transcript $transcript too short.\n";

            $headToSkip = int ( $headToSkip/2 );
            $tailToSkip = int ( $tailToSkip/2 );
            $trimed_last = $len - $tailToSkip;
            last if ( ( $headToSkip == 0 ) && ( $tailToSkip == 0 ) );
        }
        ##Basedensity###
        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );
        if ($scalingFactor > 1) {

            ###slice transcript start
            print OUT2 $transcript, "\t", $len, "\tbaseDensity\t", $rpkm, "\t" ;
            if ( $len < $windowSize ){
                $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );
                print OUT2 $scalingFactor;
                if ($scalingFactor == 0 ){
                    for ( my $idx = 1; $idx <= $len; $idx++){
                        print OUT2 "\t", sprintf ( "%.3f", $baseDensity[$idx]) ;
                    }
                }
                else{
                    for ( my $idx = 1; $idx <= $len; $idx++ )  {
                        print OUT2 "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                    }
                }
                
                print OUT2 "\n";
            }
            else{
                my @scalingFactors = ();
                my %normalizedBaseDensity = ();
                my @tmpSlicedBaseDensity = ();
                my $windowHeadToSkip = 0;
                my $windowTrimmedLast = $windowSize - 1;
                my $windowCount = 0;
                my $windowStart = 0;
                my $windowEnd = 0;
                my $tran_idx = 0;
                for ( my $i = 1; $i <= $len - $windowSize + 1; $i = $i + $stepSize ){
                    $windowStart = $i;
                    $windowEnd = $i + $windowSize - 1;
                     @tmpSlicedBaseDensity = @baseDensity[$windowStart..$windowEnd] ;# slice transcript
                     $windowCount++ ;
                     # print $transcript, "\twindow No.$windowCount\tstart at $windowStart end at $windowEnd \n";
                     my $tmpScalingFactor = 1 ;### initialize as 1 
                     
                     if ( $windowStart < $headToSkip ){ 
                         $windowHeadToSkip = $headToSkip - $windowStart ; 
                     }
                     else{
                         $windowHeadToSkip = 0;
                     } # chomp the head and tail when calculate scalingFactor
                     if ( $windowEnd > $trimed_last){ 
                         $windowTrimmedLast = $windowSize - 1 + ($trimed_last - $windowEnd); 
                     }
                     else{
                        $windowTrimmedLast = $windowSize - 1;
                     }

                    $tmpScalingFactor = &calcScalingFactor ( \@tmpSlicedBaseDensity, $windowHeadToSkip, $windowTrimmedLast, $normalizeMethod) ;
                    $tmpScalingFactor = ( $tmpScalingFactor / $scalingForm ) ;
                    # print "scalingFactor: $tmpScalingFactor using window start at $windowHeadToSkip end at $windowTrimmedLast\n";
                    
                    push @scalingFactors, $tmpScalingFactor ;

                    if ( $tmpScalingFactor == 0 ){
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, $tmpSlicedBaseDensity[$win_idx] ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [$tmpSlicedBaseDensity[$win_idx]] ; 
                            }
                            #push @normalizedBaseDensity, $tmpSlicedBaseDensity[$idx] ;
                        }
                    }
                    else{
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1 ;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, ($tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [ ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ] ; 
                            }
                            #push @normalizedBaseDensity, ($tmpSlicedBaseDensity[$idx] / $tmpScalingFactor ) ;
                        }
                    }
                    
                    # print scalar( keys %normalizedBaseDensity)."\n" ;
                }
                #final window check
                if ( $windowEnd < $len ){
                    my $leftLen = $len - $windowEnd;
                    $windowStart = $len - $windowSize + 1;
                    $windowEnd = $len;
                    @tmpSlicedBaseDensity = @baseDensity[$windowStart..$windowEnd];
                    $windowCount++;
                    # print $transcript, "\twindow No.$windowCount\tstart at $windowStart end at $windowEnd \n";
                    my $tmpScalingFactor = 1 ;
                    if ( $windowStart < $headToSkip ){ 
                        $windowHeadToSkip = $headToSkip - $windowStart ; 
                    }
                    else{
                        $windowHeadToSkip = 0;
                    } # chomp the head and tail when calculate scalingFactor
                    if ( $windowEnd > $trimed_last){ 
                        $windowTrimmedLast = $windowSize - 1 + ($trimed_last - $windowEnd); 
                    }
                    else{
                        $windowTrimmedLast = $windowSize - 1;
                    }

                    $tmpScalingFactor = &calcScalingFactor ( \@tmpSlicedBaseDensity, $windowHeadToSkip, $windowTrimmedLast, $normalizeMethod) ;
                    $tmpScalingFactor = ( $tmpScalingFactor / $scalingForm ) ;
                    # print "scalingFactor: $tmpScalingFactor using window start at $windowHeadToSkip end at $windowTrimmedLast\n";
                    push @scalingFactors, $tmpScalingFactor;
                    if ( $tmpScalingFactor == 0){
                        for ( my $win_idx = $windowSize - $leftLen; $win_idx < $windowSize; $win_idx++){
                            $tran_idx = $win_idx + $windowStart - 1 ;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, $tmpSlicedBaseDensity[$win_idx] ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [$tmpSlicedBaseDensity[$win_idx]] ; 
                            }
                            #push @normalizedBaseDensity, $tmpSlicedBaseDensity[$idx] ;
                        }
                    }
                    else{
                        for ( my $win_idx = $windowSize - $leftLen; $win_idx < $windowSize; $win_idx++){
                            $tran_idx = $win_idx + $windowStart - 1 ;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [ ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ] ; 
                            }
                            #push @normalizedBaseDensity, ( $tmpSlicedBaseDensity[$idx] / $tmpScalingFactor );
                        }
                    }
                    
                }

                die "error in length of normalized basedensity: ".scalar( keys %normalizedBaseDensity)."\n" if (scalar( keys %normalizedBaseDensity) != $len );
                print OUT2 join(",", @scalingFactors);
                my @finalBaseDensity = ();
                for (my $idx = 0; $idx < scalar( keys %normalizedBaseDensity ) ; $idx++ ){
                    my $value = eval( join( "+", @{$normalizedBaseDensity{$idx}} ) ) / scalar( @{$normalizedBaseDensity{$idx}} );   
                    push @finalBaseDensity, $value ;
                    print OUT2 "\t", sprintf ( "%.3f", $value ) ;
                }
                print OUT2 "\n";
            }
        }
=head
                if (( defined $stepSize ) and ( $stepSize > 1 ) ){
                    my $smoothedBaseDensity = &smooth( scalar(@normalizedBaseDensity), \@normalizedBaseDensity, $stepSize);
                    #print scalar(@$smoothedBaseDensity);
                    for (my $idx = 0; $idx < scalar(@$smoothedBaseDensity); $idx++){
                        print OUT2 "\t", sprintf ( "%.3f", ( $smoothedBaseDensity->[$idx] ) );
                    }
                }
                else{
                    for ( my $idx = 0; $idx <= $#normalizedBaseDensity ; $idx++){
                     print OUT2 "\t", sprintf ( "%.3f", ( $normalizedBaseDensity[$idx] ) );
                    }
                }
=cut
                
            
            ####slice transcript end
        
        ####normal one
=head
        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );

        if ( $scalingFactor > 1 ) {
            $scalingFactor = ( $scalingFactor / $scalingForm );

            print OUT $transcript, "\t", $len, "\tbaseDensity\t", $rpkm, "\t", $scalingFactor;
            if ( $normalizeMethod =~ /log/i ) {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", &log2 ( $baseDensity[$idx]/$scalingFactor + 1 ) );
                }
            }
            else {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                }
            }
            print OUT "\n";
        }
=cut
        ####RTstop####

        $line = <IN>;
        chomp $line;
        ( $transcript, $len, $rpkm, @baseDensity ) = split ( /\t/, $line );
=head
        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );

        if ( $scalingFactor > 1 ) {
            $scalingFactor = ( $scalingFactor / $scalingForm );

            print OUT $transcript, "\t", $len, "\tRTstop\t", $rpkm, "\t", $scalingFactor;
            if ( $normalizeMethod =~ /log/i ) {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", &log2 ( $baseDensity[$idx]/$scalingFactor + 1 ) );
                }
            }
            else {
                for ( my $idx = 1; $idx <= $len; $idx++ )  {
                    print OUT "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                }
            }
            print OUT "\n";
        }
=cut
        $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );
        if ($scalingFactor > 1) {
            #slice transcript
            print OUT2 $transcript, "\t", $len, "\tRTstop\t", $rpkm, "\t" ;
            if ( $len < $windowSize ){
                $scalingFactor = &calcScalingFactor ( \@baseDensity, $headToSkip, $trimed_last, $normalizeMethod );
                print OUT2 $scalingFactor;
                if ( $scalingFactor == 0 ){
                    for ( my $idx = 1; $idx <= $len; $idx++ )  {
                        print OUT2 "\t", sprintf ( "%.3f", $baseDensity[$idx] );
                    }
                }
                else{
                    for ( my $idx = 1; $idx <= $len; $idx++ )  {
                        print OUT2 "\t", sprintf ( "%.3f", ( $baseDensity[$idx]/$scalingFactor ) );
                    }
                }
                print OUT2 "\n";
            }
            else{
                my @scalingFactors = ();
                my %normalizedBaseDensity = ();
                my @tmpSlicedBaseDensity = ();
                my $windowHeadToSkip = 0;
                my $windowTrimmedLast = $windowSize - 1;
                my $windowCount = 0;
                my $windowStart = 0;
                my $windowEnd = 0;
                my $tran_idx = 0;
                for ( my $i = 1; $i <= $len - $windowSize + 1; $i = $i + $stepSize ){
                    $windowStart = $i;
                    $windowEnd = $i + $windowSize - 1;
                     @tmpSlicedBaseDensity = @baseDensity[$windowStart..$windowEnd] ;# slice transcript
                     $windowCount++ ;
                     # print $transcript, "\twindow No.$windowCount\tstart at $windowStart end at $windowEnd \n";
                     my $tmpScalingFactor = 1 ;### initialize as 1 
                     
                     if ( $windowStart < $headToSkip ){ 
                         $windowHeadToSkip = $headToSkip - $windowStart ; 
                     }
                     else{
                         $windowHeadToSkip = 0;
                     } # chomp the head and tail when calculate scalingFactor
                     if ( $windowEnd > $trimed_last){ 
                         $windowTrimmedLast = $windowSize - 1 + ($trimed_last - $windowEnd); 
                     }
                     else{
                        $windowTrimmedLast = $windowSize - 1;
                     }

                    $tmpScalingFactor = &calcScalingFactor ( \@tmpSlicedBaseDensity, $windowHeadToSkip, $windowTrimmedLast, $normalizeMethod) ;
                    $tmpScalingFactor = ( $tmpScalingFactor / $scalingForm ) ;
                    # print "scalingFactor: $tmpScalingFactor using window start at $windowHeadToSkip end at $windowTrimmedLast\n";
                    
                    push @scalingFactors, $tmpScalingFactor ;
                    if ($tmpScalingFactor == 0){
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, $tmpSlicedBaseDensity[$win_idx] ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [$tmpSlicedBaseDensity[$win_idx]] ; 
                            }
                            #push @normalizedBaseDensity, $tmpSlicedBaseDensity[$idx] ;
                        }
                    }
                    else{
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1 ;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [ ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ] ; 
                            }
                            #push @normalizedBaseDensity, ($tmpSlicedBaseDensity[$idx] / $tmpScalingFactor ) ;
                        }
                    }
                    
                    # print scalar( keys %normalizedBaseDensity)."\n" ;
                }
                #final window check
                if ( $windowEnd < $len ){
                    my $leftLen = $len - $windowEnd;
                    $windowStart = $len - $windowSize + 1;
                    $windowEnd = $len;
                    @tmpSlicedBaseDensity = @baseDensity[$windowStart..$windowEnd];
                    $windowCount++;
                    # print $transcript, "\twindow No.$windowCount\tstart at $windowStart end at $windowEnd \n";
                    my $tmpScalingFactor = 1 ;
                    if ( $windowStart < $headToSkip ){ 
                        $windowHeadToSkip = $headToSkip - $windowStart ; 
                    }
                    else{
                        $windowHeadToSkip = 0;
                    } # chomp the head and tail when calculate scalingFactor
                    if ( $windowEnd > $trimed_last){ 
                        $windowTrimmedLast = $windowSize - 1 + ($trimed_last - $windowEnd); 
                    }
                    else{
                        $windowTrimmedLast = $windowSize - 1;
                    }

                    $tmpScalingFactor = &calcScalingFactor ( \@tmpSlicedBaseDensity, $windowHeadToSkip, $windowTrimmedLast, $normalizeMethod) ;
                    $tmpScalingFactor = ( $tmpScalingFactor / $scalingForm ) ;
                    # print "scalingFactor: $tmpScalingFactor using window start at $windowHeadToSkip end at $windowTrimmedLast\n";
                    push @scalingFactors, $tmpScalingFactor;
                    if ($tmpScalingFactor == 0){
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, $tmpSlicedBaseDensity[$win_idx] ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [$tmpSlicedBaseDensity[$win_idx]] ; 
                            }
                            #push @normalizedBaseDensity, $tmpSlicedBaseDensity[$idx] ;
                        }
                    }
                    else{
                        for ( my $win_idx = 0; $win_idx < $windowSize; $win_idx++ ){
                            $tran_idx = $win_idx + $windowStart - 1 ;
                            if ( defined $normalizedBaseDensity{$tran_idx} ){
                                 push @{$normalizedBaseDensity{$tran_idx}}, ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ;
                            }
                            else{
                                $normalizedBaseDensity{$tran_idx} = [ ( $tmpSlicedBaseDensity[$win_idx] / $tmpScalingFactor ) ] ; 
                            }
                            #push @normalizedBaseDensity, ($tmpSlicedBaseDensity[$idx] / $tmpScalingFactor ) ;
                        }
                    }
                    
                }
                die "error in length of normalized basedensity: ".scalar( keys %normalizedBaseDensity)."\n" if (scalar(keys %normalizedBaseDensity) != $len);
                print OUT2 join(",", @scalingFactors);
                my @finalRTstop = ();
                for (my $idx = 0; $idx < scalar( keys %normalizedBaseDensity) ; $idx++ ){
                    my $value = eval( join( "+", @{$normalizedBaseDensity{$idx}} ) ) / scalar( @{$normalizedBaseDensity{$idx}} );   
                    push @finalRTstop, $value ;
                    print OUT2 "\t", sprintf ( "%.3f", $value ) ;
                }
                print OUT2 "\n";
            }
        }
=head
                if (( defined $stepSize ) and ($stepSize > 1)){
                    my $smoothedBaseDensity = &smooth( scalar(@normalizedBaseDensity), \@normalizedBaseDensity, $stepSize);
                    for (my $idx = 0; $idx < scalar(@$smoothedBaseDensity) ; $idx++){
                        print OUT2 "\t", sprintf ( "%.3f", ( $smoothedBaseDensity->[$idx] ) );
                    }
                }
                else{
                    for ( my $idx = 0; $idx <= $#normalizedBaseDensity ; $idx++){
                     print OUT2 "\t", sprintf ( "%.3f", ( $normalizedBaseDensity[$idx] ) );
                    }
                }
=cut
                
            
    }
    close IN;
    #close OUT;
    close OUT2;
    print STDERR "Normalize base density from file $baseDensityFile successfully, $transCount transcripts processed\n\t", `date`;
    
    1;
}

sub calcScalingFactor  {
    my ( $ref_array, $start, $end, $normalizeMethod ) = @_;

    my @sortedTrimed_array = sort { $b <=> $a } ( @{$ref_array}[$start..$end] );

    my $len = scalar(@sortedTrimed_array); 
    my @rangeOfSelection = ();
    if ( $normalizeMethod =~ /smart/i ) {
        ## smartly skip those non-zero elements
        for ( my $idx = 0; $idx < $len; $idx++ )  {
            last if ( $sortedTrimed_array[$idx] <= 0 );
            push @rangeOfSelection, $sortedTrimed_array[$idx];
        }
    }
    else {
        my $selectStart = 0; 
        my $selectEnd = $len - 1;
        if ( $normalizeMethod =~ /upper/i ) {
            ## get the upper half
            $selectEnd = int ( $len / 2 ) - 1;
        }
        elsif ( $normalizeMethod =~ /quartile(\d+)/i ) {
            ## get a specified quarter
            $selectStart = int ( ($1-1)*$len / 4 );
            $selectEnd = int ( $1*$len / 4 ) - 1;
        }
        elsif ( $normalizeMethod =~ /decile(\d+)/i ) {
            ## get a specified decile
            $selectStart = int ( ($1-1)*$len / 10 );
            $selectEnd = int ( $1*$len / 10 ) - 1;
        }
        elsif ( $normalizeMethod =~ /vigintile(\d+)/i ) {
            ## get a specified decile
            $selectStart = int ( ($1-1)*$len / 20 );
            $selectEnd = int ( $1*$len / 20 ) - 1;
        }

        @rangeOfSelection = @sortedTrimed_array[$selectStart..$selectEnd];
    }
    $len = scalar(@rangeOfSelection); 

    my $scalingFactor = 1;
    if ( $normalizeMethod =~ /median/i ) {
        my $median = 1;
        if ( $len % 2 == 0 ) {  $median = ($rangeOfSelection[$len/2-1] + $rangeOfSelection[$len/2]) /2;  }
        else {  $median = $rangeOfSelection[($len-1)/2];  }
        $scalingFactor = $median;
    }
    elsif ( $normalizeMethod =~ /mean/i ) {
        my $mean = eval ( join "+", @rangeOfSelection );
        $mean /= scalar (@rangeOfSelection);
        $scalingFactor = $mean;
    }
    elsif ( $normalizeMethod =~ /peak/i ) {
        $scalingFactor = $rangeOfSelection[0];
    }

    return $scalingFactor;
}

sub log2 {
    my $num = shift;
    return log($num)/log(2);
}
=head
sub smooth {
    my $len = shift;
    my $ref_array = shift;
    my $smoothWindowSize = shift;
    my @smoothedArray = ();
    my $halfWindowSize = int( $smoothWindowSize / 2 );
    my @tmpSmoothWindow = ();
    ##the first halfWindowSize bases
    my $start = 0;
    my $end = 0;
    for (my $i = 0 ; $i < $halfWindowSize; $i++){
        $start = 0;
        $end = $i + $halfWindowSize;
        @tmpSmoothWindow = @{$ref_array}[$start..$end];
        @smoothedArray[$i] = eval( join( "+", @tmpSmoothWindow)) / scalar( @tmpSmoothWindow );

    }
    #middle
    
    for (my $i = $halfWindowSize; $i < $len - $halfWindowSize; $i++){
        $start = $i - $halfWindowSize;
        $end = $i + $halfWindowSize;
        @tmpSmoothWindow = @{$ref_array}[$start..$end];
        @smoothedArray[$i] = eval( join( "+", @tmpSmoothWindow)) / scalar( @tmpSmoothWindow );
    } 
    #the last smoothWindowSize bases
    for (my $i = $len - $halfWindowSize; $i < $len; $i++){
        $start = $i - $halfWindowSize;
        $end = $len - 1 ;
        @tmpSmoothWindow = @{$ref_array}[$start..$end];
        @smoothedArray[$i] = eval( join( "+", @tmpSmoothWindow)) / scalar( @tmpSmoothWindow );
    }   
    return \@smoothedArray;
    1;
}
=cut


