#!/usr/bin/env perl

# job_make_rnaseq_renaming_25may2018.pl -- Erich Schwarz <ems394@cornell.edu>, 5/25/2018.
# Purpose: take data plain-text (from Excel) as input; generate reliable renaming script for RNA-seq read files.

use strict;
use warnings;
use autodie;

use File::Basename;

my $input_file_list = $ARGV[0];
my $data_table      = $ARGV[1];

my $data_ref;

my %seen = ();

open my $INLIST, '<', $input_file_list;
while (my $input = <$INLIST>) {
    chomp $input;
    my $basename = basename($input);
    if ( exists $data_ref->{'basename'}->{$basename}->{'fullname'} ) {
        die "Redundant full file names for $basename: $input and $data_ref->{'basename'}->{$basename}->{'fullname'}\n";
    }
    $data_ref->{'basename'}->{$basename}->{'fullname'} = $input;
}
close $INLIST;

open my $DATA, '<', $data_table;
while (my $input = <$DATA>) {
    chomp $input;

    # Sample input line:
    # T727    TWF154  T727-10hr       10 hr after N2 added    LNM     1       32-T727-10-1_S5_R1_001.fastq.gz N414-20170619-TruSeq_RNA

    if ( $input =~ /\A T727 \t [^\t]* \t (\S+) \t [^\t]* \t [^\t]* \t (\d+) \t (\S+) \t [^\t]* \z/xms ) { 
        my $condition     = $1;
        my $replicate     = $2;
        my $orig_filename = $3;
        my $orig_basename = basename($orig_filename);
        my $read_no       = 0;

        # Sample basename:
        # 28-T727-24-3_S9_R2_001.fastq.gz

        if ( $orig_basename =~ /\A \S+ _ (R[1-2]) _ \d+ \.fastq\.gz \z/xms ) {
            $read_no = $1;
        }
        else { 
            die "Cannot identify read number in basename: $orig_basename\n";
        }
        if (! exists $data_ref->{'basename'}->{$orig_basename}->{'fullname'} ) {
            die "Cannot identify full file name for: $orig_basename\n";
        }

        my $new_filename  = "$condition" . ".rep_" . "$replicate" . "_" . "$read_no.raw.fastq";
        if (-e $new_filename) {
            die "$new_filename already exists in this directory\n";
        }

        if (! $seen{$condition} ) {
            print "\n";
            $seen{$condition} = 1;
        }
        if (! $seen{$orig_basename} ) {
            if ( $seen{$new_filename} ) {
                print "WARNING: ";
            }
            else {
                $seen{$new_filename} = 1;
            }
            my $full_orig_filename = $data_ref->{'basename'}->{$orig_basename}->{'fullname'};
            print "zcat $full_orig_filename > $new_filename ;\n";
            $seen{$orig_basename} = 1;
        }
    }
    elsif ( $input !~ /\A Strain \t /xms ) { 
        die "Cannot parse: $input\n";
    }
}
print "\n";
close $DATA ;


