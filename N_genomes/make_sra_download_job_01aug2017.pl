#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $srr    = q{};
my $prefix = q{};

my $header = 
    "#!/bin/bash -login\n"
    . "#PBS -l walltime=048:00:00\n"
    . "#PBS -l nodes=1:ppn=1\n"
    . "#PBS -l mem=32gb\n"
    . "#PBS -N GENERIC_COMMAND_NAME.sh\n"
    . "#PBS -q main\n"
    . '#PBS -M ems394@cornell.edu' . "\n"
    . "#PBS -m abe\n"
    . "#PBS -A ged\n"
    . "#PBS -r n\n"
    . "#PBS -V\n"
    ;

while (my $infile = <>) {
    chomp $infile;

    if (! -r $infile) { 
        die "Cannot read input file: $infile\n";
    }

    my $dirname  = dirname($infile);

    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;

        if ( $input =~ /\A .+ (SRR \d+) \s* \z/xms ) { 
            $srr = $1;
        }
        else { 
            die "Can't parse input line: $input\n";
        }
        if ( $srr =~ /\A (SRR\d{3}) \d+ \z/xms ) { 
            $prefix = $1;
        }
        else {
            die "Can't parse putative SRR accession number : $srr\n";
        }

        print "$header" if $header;
        $header = q{};

        print "cd $dirname ;\n";
        print 'wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/' . $prefix . q{/} . $srr . q{/} . "$srr.sra ;\n";
        print "fastq-dump --split-3 $srr.sra ;\n";
    }
    close $INFILE;
}

