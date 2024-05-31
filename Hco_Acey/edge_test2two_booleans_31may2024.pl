#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $infile_list = q{};

$infile_list = $ARGV[0] if $ARGV[0];

if (! $infile_list ) {
    die "Format: edge_test2boolean_31may2024.pl [list of input edgeR files] => [2 autonamed Boolean annotation file -- upreg. vs. downreg.]\n"
}

open my $INFILE_LIST, '<', $infile_list;
while (my $infile = <$INFILE_LIST>) {
    chomp $infile;
    if (! -r $infile ) {
        die "Cannot read infile: $infile\n";
    }

    my $basename = basename($infile);
    my $stem     = q{};
    if ( $basename =~ /\A (\S+)_edgeR_exactTest_\d+\.\d+\.\d+\.\d+\.tsv.txt \z/xms ) {
        $stem = $1;
    }
    else {
        die "Cannot parse basename: $basename\n";
    }

    my $up_outfile = "$stem.upreg_annot.tsv.txt";
    $up_outfile    = safename($up_outfile);

    my $down_outfile = "$stem.downreg_annot.tsv.txt";
    $down_outfile    = safename($down_outfile);

    open my $UP_OUTFILE,   '>', $up_outfile;
    open my $DOWN_OUTFILE, '>', $down_outfile;

    open my $INFILE, '<', $infile;
    my $up_header = "Gene\t$stem.upreg";
    my $down_header = "Gene\t$stem.downreg";
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input !~ /\A (?: \S+ \t){4} \S+ \z/xms ) {
            die "Mis-tabbed line in input file $infile: $input\n";
        }
        elsif ( $input =~ /\A Gene \t (\S+)\.logFC \t (\S+)\.logCPM \t (\S+)\.PValue \t (\S+)\.FDR \z/xms ) {
            my $stem1 = $1;
            my $stem2 = $2;
            my $stem3 = $3;
            my $stem4 = $4;
            if ( ( $stem1 ne $stem ) or ( $stem2 ne $stem ) or ( $stem3 ne $stem ) or ( $stem4 ne $stem ) ) {
                die "Mis-stemmed headers in input file $infile: $input\n";
            }
        }
        elsif ( $input =~ /\A (\S+) \t (\S+) \t \S+ \t \S+ \t (\S+) \z/xms  ) {
            my $gene  = $1;
            my $logFC = $2;
            my $fdr   = $3;

            print $UP_OUTFILE "$up_header\n" if $up_header;
            print $DOWN_OUTFILE "$down_header\n" if $down_header;
            $up_header   = q{};
            $down_header = q{};

            if ( $fdr <= 0.01) {
                if ( $logFC >= 1 ) {
                    print $UP_OUTFILE "$gene\t$stem.upreg\n";
                }
                elsif ( $logFC <= -1 ) {
                    print $DOWN_OUTFILE "$gene\t$stem.downreg\n";
                }
            }
        }
        else {
            die "From input file $infile, cannot parse: $input\n";
        }
    }
    close $INFILE;
    close $UP_OUTFILE;
    close $DOWN_OUTFILE;
}
close $INFILE_LIST;

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

