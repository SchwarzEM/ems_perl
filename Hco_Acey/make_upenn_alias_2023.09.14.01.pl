#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir
use File::Basename;

my $data_ref;

my $start_dir = getcwd;

my $accessions = q{};
my $prev_accs  = q{};
my $renaming   = q{};
my $readlist   = q{};

$accessions = $ARGV[0] if $ARGV[0];
$prev_accs  = $ARGV[1] if $ARGV[1];
$renaming   = $ARGV[2] if $ARGV[2];
$readlist   = $ARGV[3] if $ARGV[3];

if ( (! $accessions ) or (! $prev_accs ) or (! $renaming ) or (! $readlist ) ) {
    die "Format: make_upenn_alias_2023.09.02.01.pl",
        " [accession/condition table]",
        " [prev. acc./cond. table]",
        " [renaming table]",
        " [read file list]",
        " > [line-commands]",
        "\n",
        ;
}

open my $ACCESSIONS, '<', $accessions;
while (my $input = <$ACCESSIONS> ) {
    chomp $input;

    # Sample input line:
    # SAMN37203688    WT_G7

    if ( $input =~ /\A (SAMN\d+) \t (\S+) \z/xms ) {
        my $accession = $1;
        my $condition = $2;

        if ( exists $data_ref->{'accession'}->{$accession}->{'condition'} ) {
            die "Redundant mapping of accession $accession to conditions $data_ref->{'accession'}->{$accession}->{'condition'} versus $condition\n";
        }
        else {
            $data_ref->{'accession'}->{$accession}->{'condition'} = $condition;
        }
    }
    else {
        die "From accessions $accessions, cannot parse: $input\n";
    }
}
close $ACCESSIONS;

open my $PREV_ACCS, '<', $prev_accs;
while (my $input = <$PREV_ACCS> ) {
    chomp $input;

    # Sample input line:
    # SAMN37203688    WT_G7

    if ( $input =~ /\A (SAMN\d+) \t (\S+) \z/xms ) {
        my $accession = $1;
        my $prev_cond = $2;

        if ( exists $data_ref->{'prev_cond'}->{$prev_cond}->{'accession'} ) {
            die "Redundant mapping of condition $prev_cond to accessions $data_ref->{'prev_cond'}->{$prev_cond}->{'accession'} versus $accession\n";
        }
        else {
            $data_ref->{'prev_cond'}->{$prev_cond}->{'accession'} = $accession;
        }
    }
    else {
        die "From prev. accs. $prev_accs, cannot parse: $input\n";
    }
}
close $PREV_ACCS;

open my $RENAMING, '<', $renaming;
while (my $input = <$RENAMING> ) {
    chomp $input;

    # Sample input line:
    # STAT6_G14   STAT6.KO_G14
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $old_cond = $1;
        my $new_cond = $2;

        if ( exists $data_ref->{'old_cond'}->{$old_cond}->{'new_cond'} ) {
            die "Redundant mapping of old condition $old_cond",
                " to new conditions $data_ref->{'old_cond'}->{$old_cond}->{'new_cond'} versus $new_cond",
                "\n",
                ;
        }
        else {
            $data_ref->{'old_cond'}->{$old_cond}->{'new_cond'} = $new_cond;
        }
    }

    else {
        die "From renaming $renaming, cannot parse: $input\n";
    }
}
close $RENAMING;

open my $READLIST, '<', $readlist;
while (my $full_readfile = <$READLIST> ) {
    chomp $full_readfile;

    if (! -e $full_readfile ) {
        die "Apparently, this input file does not exist: $full_readfile\n";
    }

    # Sample input full readfile:
    # /ocean/projects/mcb190015p/fergusoa/Herbert_Lab_Nb_RNAseq_data/AF37_STAT6G6_WTG7_fastq/01_WT_G7_trm.fastq.gz
    if ( $full_readfile =~ /\A \/ \S+ \/ \d+ _ (\S+) _trm\.fastq\.gz /xms ) { 
        my $prev_cond  = $1;

        # This is completely arbitrary, but it's necessary to deal with an inconsistently named replicate set
        #    that is now locked into an official BioProject name.
        if ( $prev_cond eq 'STAT6_G6' ) {
            $prev_cond = 'STAT6KO_G6';
        }

        if (! exists $data_ref->{'prev_cond'}->{$prev_cond}->{'accession'} ) {
            die "Cannot map previous condition $prev_cond to accession\n";
        }
        my $accession = $data_ref->{'prev_cond'}->{$prev_cond}->{'accession'};
        my $condition = $data_ref->{'accession'}->{$accession}->{'condition'} ;
        my $sub_dir = "$accession.$condition". '_02';

        my $work_dir = catdir($start_dir, $sub_dir);

        my $outfile = basename($full_readfile);
        if ( $outfile =~ /\A (\d+ _) (\S+) (_trm\.fastq\.gz) \z/xms ) {
            my $prefix = $1;
            my $tag    = $2;
            my $suffix = $3;
            if (! exists $data_ref->{'old_cond'}->{$tag}->{'new_cond'} ) {
                die "Cannot rename tag: $tag\n";
            }
            else {
                $tag = $data_ref->{'old_cond'}->{$tag}->{'new_cond'};
            }
            $outfile = "$prefix$tag$suffix";
        }
        else {
            die "Cannot rename intended outfile: $outfile\n";
        }

        print "mkdir $work_dir ;\n";
        print "cd $work_dir ;\n";
        print "ln -s $full_readfile $outfile ;\n";
    }
    else {
        die "Cannot parse input: $full_readfile\n";
    }
}
close $READLIST;
