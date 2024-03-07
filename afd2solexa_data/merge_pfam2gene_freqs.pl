#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

if (! @ARGV ) {
    die "Format: merge_pfam2gene_freqs.pl [one or more pfam2gene_freq infiles] > [one merged pfam2gene_freq outfile]\n";
}

my $data_ref;

my @infiles     = @ARGV;
my $infile_text = join "\t", @infiles;
my $header      = "Motif\t$infile_text";

if (@infiles) {
    foreach my $infile (@infiles) {
        open my $INFILE, '<', $infile;
        while ( my $input = <$INFILE> ) {
            chomp $input;
            if ( $input =~ /\A ([^\t]+) \t (\d+) \z/xms ) {
                my $motif = $1;
                my $count = $2;
                $data_ref->{'motif'}->{$motif}->{'infile'}->{$infile}->{'count'} = $count;
            }
        }
        close $INFILE;
    }
}

my @motifs = sort keys %{ $data_ref->{'motif'} };

foreach my $motif (@motifs) {
    my $output = $motif;
    foreach my $infile (@infiles) {
        my $count = 0;
        if ( exists $data_ref->{'motif'}->{$motif}->{'infile'}->{$infile}->{'count'} ) {
            $count = $data_ref->{'motif'}->{$motif}->{'infile'}->{$infile}->{'count'};
        }
        $output = "$output\t$count";
    }
    print "$header\n" if $header;
    $header = q{};
    print "$output\n";
}
