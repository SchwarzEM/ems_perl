#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use File::Spec::Functions;

my $list = $ARGV[0];

open my $LIST, '<', $list;
while (my $input = <$LIST>) {
    chomp $input;
    if (! -r $input) {
        die "From input list file $list, cannot read putative file: $input\n";
    }
    my $basename = basename($input);
    $basename    =~ s/\.txt\z//;

    my $script   = 'job_meme_' . $basename . '_2016.03.07.sh';
    $script      = safename($script);

    my $full_script = catfile('/mnt/home/emsch/work/2015/adrienne/ncDNA_motifs', $script);

    # Every time this runs, it starts a fresh work directory, making it easier to do re-runs.
    my $work_dir = $basename . '_dir';
    $work_dir    = catfile('/mnt/home/emsch/work/2015/adrienne/ncDNA_motifs', $basename);
    $work_dir    = safename($work_dir);

    my $promseqs  = $basename . '_500trans_dna.fa';
    my $memedir   = $basename . '_500trans_meme_2016.03.07';
    my $tomtomdir = $basename . '_500trans_tomtom_2016.03.07';
    my $meme_file = catfile($memedir, 'meme.txt');

    open my $SCRIPT, '>', $script;

    print $SCRIPT '#!/bin/bash -login', "\n";
    print $SCRIPT '#PBS -l walltime=004:00:00', "\n";
    print $SCRIPT '#PBS -l nodes=1:ppn=8', "\n";
    print $SCRIPT '#PBS -l mem=16gb', "\n";
    print $SCRIPT '#PBS -N ', "$script\n";
    print $SCRIPT '#PBS -q main', "\n";
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
    print $SCRIPT '#PBS -m abe', "\n";
    print $SCRIPT '#PBS -A ged', "\n";
    print $SCRIPT '#PBS -r n', "\n";
    print $SCRIPT '#PBS -V', "\n";
    print $SCRIPT "cd /mnt/home/emsch/work/2015/adrienne/ncDNA_motifs ;\n";
    print $SCRIPT "module load MEME/4.11.1 ;\n";
    print $SCRIPT "mkdir $work_dir ;\n";
    print $SCRIPT "cd $work_dir ;\n";

    print $SCRIPT "extract_fasta_subset.pl -w -l $input ";
    print $SCRIPT "-f /mnt/home/emsch/work/2015/adrienne/ncDNA_motifs/TAIR10_upstream_500_translation_start_20101028.fa 1>";
    print $SCRIPT "$promseqs ;\n";

    print $SCRIPT "meme $promseqs -evt 0.1 -o $memedir -dna -mod zoops -nmotifs 20 -minw 6 -maxw 12 ";
    print $SCRIPT "-bfile /mnt/home/emsch/work/2015/adrienne/ncDNA_motifs/TAIR10_500nt_trans_markov1 -revcomp -maxsize 10000000 -p 8 ;\n";

    print $SCRIPT "tomtom -o $tomtomdir -thresh 0.1 ";
    print $SCRIPT "-bfile /mnt/home/emsch/work/2015/adrienne/ncDNA_motifs/TAIR10_500nt_trans_markov1 ";
    print $SCRIPT "$meme_file ";
    print $SCRIPT "/mnt/home/emsch/meme_4.11.1/db/{ARABD,CIS-BP,EUKARYOTE,FLY,HUMAN,JASPAR,MALARIA,MOUSE,TFBSshape,WORM,YEAST}/*.meme ;\n";

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


