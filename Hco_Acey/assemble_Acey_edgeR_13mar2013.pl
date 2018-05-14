#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my @input_files = @ARGV;

foreach my $input_file (@input_files) { 
    my $input_file_basename = basename $input_file;
    if ( $input_file_basename =~ /\A Acey_counts\. (\S+ \.vs\. \S+) \.DGEList\.exactTest\.topTags\.txt \z/xms ) {  
        my $stem        = $1;
        $stem           = $stem . '.edgeR';
        my $output_file = $stem . '.txt';
        $output_file    = safename($output_file);
        my $header      = "Gene\t$stem.logFC\t$stem.q-value\n";

        open my $OUTPUT_FILE, '>', $output_file or die "Can't open output file $output_file: $!";
        open my $INPUT_FILE, '<', $input_file or die "Can't open input file $input_file: $!";

        while ( my $input = <$INPUT_FILE> ) { 
            chomp $input;
            # Sample input lines:
            # 
            # [\t] "logFC" "logCPM" "PValue" "FDR"
            # "Acey_2012.08.05_0038.g3635" 19.8178817085073 13.4305837864642 4.34698549465812e-19 9.51011666162684e-15
            if ( $input =~ /\A \" (\S+) \" \s+ (\S+) \s+ \S+ \s+ \S+ \s+ (\S+) \z/xms ) {
                my $gene  = $1;
                my $logFC = $2;
                my $q_val = $3;
                if ( ( looks_like_number($logFC) ) and ( looks_like_number($q_val) ) ) { 
                    print $OUTPUT_FILE $header if $header;
                    $header = q{};
                    print $OUTPUT_FILE "$gene\t$logFC\t$q_val\n";
                }
            }
        }

        close $INPUT_FILE or die "Can't close filehandle to input file $input_file: $!";
        close $OUTPUT_FILE or die "Can't close filehandle to output file $output_file: $!";
    }
    else { 
        die "Can't parse input file $input_file\n";
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

