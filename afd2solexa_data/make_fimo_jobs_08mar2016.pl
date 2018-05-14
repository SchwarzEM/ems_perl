#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;

my $list = $ARGV[0];

open my $LIST, '<', $list;
while (my $input = <$LIST>) {
    chomp $input;

    # Enforce reality of the target MEME file
    if (! -r $input) {
        die "From input list file $list, cannot read putative file: $input\n";
    }

    # Enforce likely usability of the target MEME file:
    if ( $input !~ / \/meme\.txt \z/xms ) { 
        die "From input list file $list, input file is unlikely to work for FIMO: $input\n";
    }

    my $work_dir = q{};
    # Do *not* safename $work_dir because it is supposed to already exist, and we are going to go back into it!
    my $stem     = q{};

    if ( $input =~ /\A ( \/mnt\/home\/emsch\/work\/2015\/adrienne\/ncDNA_motifs\/ ([^\/\s]+) ) \/ /xms ) { 
        $work_dir = $1;
        $stem     = $2;
    }
    else {
        die "Cannot parse input file name: $input\n";
    }

    my $script   = 'job_fimo_' . $stem . '_2016.03.08.sh';
    $script      = safename($script);

    my $full_script = catfile('/mnt/home/emsch/work/2015/adrienne/ncDNA_motifs', $script);

    my $fimo_dir  = $stem . '_500trans_fimo_2016.03.07';
    $fimo_dir     = safename($fimo_dir);

    # strictly speaking unnecessary, but it makes it easier to keep things straight, so do this:
    my $meme_file = $input;

    open my $SCRIPT, '>', $script;

    print $SCRIPT '#!/bin/bash -login', "\n";
    print $SCRIPT '#PBS -l walltime=001:00:00', "\n";
    print $SCRIPT '#PBS -l nodes=1:ppn=1', "\n";
    print $SCRIPT '#PBS -l mem=16gb', "\n";
    print $SCRIPT '#PBS -N ', "$script\n";
    print $SCRIPT '#PBS -q main', "\n";
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
    print $SCRIPT '#PBS -m abe', "\n";
    print $SCRIPT '#PBS -A ged', "\n";
    print $SCRIPT '#PBS -r n', "\n";
    print $SCRIPT '#PBS -V', "\n";
    print $SCRIPT "cd $work_dir ;\n";
    print $SCRIPT "module load MEME/4.11.1 ;\n";
    print $SCRIPT "fimo --bgfile /mnt/home/emsch/work/2015/adrienne/ncDNA_motifs/TAIR10_500nt_trans_markov1 ";
    print $SCRIPT "--max-stored-scores 10000000 -o $fimo_dir --thresh 1e-4 $meme_file ";
    print $SCRIPT "/mnt/home/emsch/work/2015/adrienne/ncDNA_motifs/TAIR10_upstream_500_translation_start_20101028.contigs.fa ;\n";
    print $SCRIPT "echo \"Finished run of $full_script.\" > $full_script.log ;\n";

    close $SCRIPT;
}

close $LIST;


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


