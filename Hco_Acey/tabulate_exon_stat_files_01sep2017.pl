#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $data_ref;

my @infiles = @ARGV;

my $header = "Data_source"
             . "\tTx_count"
             . "\tExon_count"
             . "\tIntron_count"
             . "\tMean_exon_sum_len"
             . "\tMean_intron_sum_len"
             . "\tMedian_exon_sum_len"
             . "\tMedian_intron_sum_len"
             . "\tMean_exon_len"
             . "\tMean_intron_len"
             . "\tMedian_exon_len"
             . "\tMedian_intron_len"
             ;

foreach my $infile (@infiles) {
    open my $INFILE, '<', $infile;

    my $basename = basename($infile);
    $basename =~ s/_max_iso_stats\.txt\z//;

    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input !~ /\A Tx_count \t/xms ) { 
            $input = "$basename\t$input";
            $data_ref->{'input_source'}->{$basename}->{'data'} = $input;
        }
    }
}

my @sources = sort keys %{ $data_ref->{'input_source'} };

foreach my $source (@sources) {
    my $data = $data_ref->{'input_source'}->{$source}->{'data'};

    print "$header\n" if $header;
    $header = q{};

    print "$data\n";
}


