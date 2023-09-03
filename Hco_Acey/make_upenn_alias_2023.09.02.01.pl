#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir

my $data_ref;

my $start_dir = getcwd;

my $accessions = q{};
my $readlist   = q{};

$accessions = $ARGV[0] if $ARGV[0];
$readlist   = $ARGV[1] if $ARGV[1];

if ( (! $accessions ) or (! $readlist ) ) {
    die "Format: make_upenn_alias_2023.09.02.01.pl [accession/condition table] [read file list] > [line-commands]\n";
}

open my $ACCESSIONS, '<', $accessions;
while (my $input = <$ACCESSIONS> ) {
    chomp $input;

    # Sample input line:
    # SAMN37203688    WT_G7

    if ( $input =~ /\A (SAMN\d+) \t (\S+) \z/xms ) {
        my $accession = $1;
        my $condition = $2;

        if ( exists $data_ref->{'condition'}->{$condition}->{'accession'} ) {
            die "Redundant mapping of condition $condition to accessions $data_ref->{'condition'}->{$condition}->{'accession'} versus $accession\n";
        }
        $data_ref->{'condition'}->{$condition}->{'accession'} = $accession;
    }
}
close $ACCESSIONS;

open my $READLIST, '<', $readlist;
while (my $input = <$READLIST> ) {
    chomp $input;

    if (! -e $input ) {
        die "Apparently, this input file does not exist: $input\n";
    }

    # Sample input:
    # /ocean/projects/mcb190015p/fergusoa/Herbert_Lab_Nb_RNAseq_data/AF37_STAT6G6_WTG7_fastq/01_WT_G7_trm.fastq.gz
    if ( $input =~ /\A \/ \S+ \/ \d+ _ (\S+) _trm\.fastq\.gz /xms ) { 
        my $condition  = $1;

        # This is completely arbitrary, but it's necessary to deal with an inconsistently named replicate set
        #    that is now locked into an official BioProject name.
        if ( $condition eq 'STAT6_G6' ) {
            $condition = 'STAT6KO_G6';
        }

        if (! exists $data_ref->{'condition'}->{$condition}->{'accession'} ) {
            die "Cannot map condition $condition to accession\n";
        }
        my $accession = $data_ref->{'condition'}->{$condition}->{'accession'};
        my $sub_dir = "$accession.$condition";

        my $work_dir = catdir($start_dir, $sub_dir);

        print "mkdir $work_dir ;\n";
        print "cd $work_dir ;\n";
        print "ln -s $input ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
close $READLIST;
