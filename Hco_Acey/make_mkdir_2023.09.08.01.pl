#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $biosamples = q{};
my $replicates = q{};
my $r1_reads   = q{};

$biosamples = $ARGV[0] if $ARGV[0];
$replicates = $ARGV[1] if $ARGV[1];
$r1_reads   = $ARGV[2] if $ARGV[2];

my $start_dir = getcwd;

my $data_ref;

if ( (! $biosamples ) or (! $replicates ) or (! $r1_reads ) ) {
    die "Format: make_mkdir_2023.09.02.01.pl [biosample table] [replicate table] [r1 read file list] > [line commands for mkdir/ln -s]\n";
}

open my $BIOSAMPLES, '<', $biosamples;
while ( my $input = <$BIOSAMPLES> ) {
    chomp $input;

    # Sample input line:
    # SAMN37187589    IJ

    if ( $input =~ /\A (SAMN\d+) \t (\S+) \z/xms ) {
        my $accession = $1;
        my $condition = $2;

        if ( exists $data_ref->{'condition'}->{$condition}->{'accession'} ) {
            die "Redundant accessions for condition $condition\n";
        }

        $data_ref->{'condition'}->{$condition}->{'accession'} = $accession;
    }
    else {
        die "From biosamples file $biosamples, cannot parse: $input\n";
    }
}
close $BIOSAMPLES;

open my $REPLICATES, '<', $replicates;
while ( my $input = <$REPLICATES> ) {
    chomp $input;

    # Sample input lines:
    # IJ      IJ_rep01        N001
    # IJ      IJ_rep02        N002
    # IJ      IJ_rep03        N003

    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
        my $condition = $1;
        my $replicate = $2;
        my $subdir    = $3;

        if ( ( exists $data_ref->{'subdir'}->{$subdir}->{'replicate'} ) or ( exists $data_ref->{'subdir'}->{$subdir}->{'condition'} ) ) {
            die "Redundant replicate and/or condition annotations for subdirectory $subdir\n";
        }

        if (! exists $data_ref->{'condition'}->{$condition}->{'accession'} ) {
            die "No accession recorded for condition $condition\n";
        }

        $data_ref->{'subdir'}->{$subdir}->{'replicate'} = $replicate;
        $data_ref->{'subdir'}->{$subdir}->{'condition'} = $condition;
    }
    else {
        die "From replicates file $replicates, cannot parse: $input\n";
    }
}
close $REPLICATES;

open my $R1_READS, '<', $r1_reads;

LOOP: while ( my $infile_1 = <$R1_READS> ) {
    chomp $infile_1;
    if ( $infile_1 =~ /\A (\S+) \/ (N\d+) \/ ( [^\s\/]+ ) \.fastq \.gz \z/xms ) {
        my $dir1_stem  = $1;
        my $subdir     = $2;
        my $file1_stem = $3;

        if ( $file1_stem !~ /R1/xms ) {
            die "File1 stem $file1_stem has no R1\n";
        }

        if ( (! exists $data_ref->{'subdir'}->{$subdir}->{'replicate'} ) or (! exists $data_ref->{'subdir'}->{$subdir}->{'condition'} ) ) {
            warn "No replicate and/or condition for subdirectory $subdir; therefore, skipping it\n";
            next LOOP;
        }

        my $replicate = $data_ref->{'subdir'}->{$subdir}->{'replicate'};
        my $condition = $data_ref->{'subdir'}->{$subdir}->{'condition'};

        if (! exists $data_ref->{'condition'}->{$condition}->{'accession'} ) {
            die "Cannot map condition $condition to accession\n";
        }
        my $accession = $data_ref->{'condition'}->{$condition}->{'accession'};

        my $infile_2   = $infile_1;
        $infile_2      =~ s/R1([^\s\/]+\.fastq\.gz)\z/R2$1/;

        my $file2_stem = $file1_stem;
        $file2_stem    =~ s/R1/R2/g;

        if ( $infile_1 eq $infile_2 ) {
            die "Failure to distinguish $infile_1 from $infile_2\n";
        }
        if ( $file1_stem eq $file2_stem ) {
            die "Failure to distinguish $file1_stem from $file2_stem\n";
        }

        if (! -r $infile_1) {
            die "Cannot find or read infile_1: $infile_1\n";
        }
        if (! -r $infile_2) {
            die "Cannot find or read infile_2: $infile_2\n";
        }

        my $outdir = catdir($start_dir, "$accession.$condition");

        # Add this to allow clearly distinct source dirs, to fix a previous sabotaged upload run.
        $outdir = $outdir . '_02';

        print "mkdir $outdir ;\n";
        print "cd $outdir ;\n";
        print "ln -s $infile_1 $replicate.R1.fastq.gz ;\n";
       	print "ln -s $infile_2 $replicate.R2.fastq.gz ;\n";
        print "\n";
    }
    else {
        die "From r1_reads file $r1_reads, cannot parse infile_1 $infile_1\n";
    }
}
close $R1_READS;

