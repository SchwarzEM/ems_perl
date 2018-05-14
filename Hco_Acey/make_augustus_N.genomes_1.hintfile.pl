#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $identifier     = q{};
my $parameters     = q{};
my @genomic_files  = ();
my $cDNA_hint_file = q{};
my $master_script  = q{};
my $aug_script     = q{};
my $utrs;

my $help;

GetOptions ( 'identifier=s' => \$identifier,
             'parameters=s' => \$parameters,
             'genome=s{,}'  => \@genomic_files,
             'cdna=s'       => \$cDNA_hint_file,
             'master=s'     => \$master_script,
             'utrs'         => \$utrs,
             'help'         => \$help, );

if ( $help or (! @genomic_files) or (! $cDNA_hint_file) or (! $master_script) ) {
    die "\n",
        "Format: make_augustus_N.genomes_1.hintfile.pl\n",
        "        --identifier|-i  [identifier to put into script names, e.g., '07aug2012']\n",
        "        --parameters|-p  [species-specific parameters to have AUGUSTUS invoke]\n",
        "        --genome|-g      [name(s) of one or more existing genomic DNA files]\n",
        "        --cdna|-c        [name existing cDNA-based GFF hint files]\n",
        "        --master|-m      [name of intended master script file]\n",
        "        --utrs|-u        [try to predict UTRs; should only do this if parameters support them]\n",
        "        --help|-h        [print this message]\n",
        "\n",
        ;
}

$master_script = safename($master_script);

open my $MASTER_SCRIPT, '>', $master_script or die "Can't open master script $master_script: $!";

print $MASTER_SCRIPT '#!/bin/bash', "\n\n", ;

foreach my $genomic_file (@genomic_files) { 
    my $gen_base  = basename $genomic_file;
    my $cDNA_base = basename $cDNA_hint_file;

    my $aug_script = 'job_augustus_' . $gen_base . '_w_' . $cDNA_base . '_hints.' . "$parameters.$identifier.sh";
    $aug_script = safename($aug_script);

    my $aug_output = $gen_base . '_w_' . $cDNA_base . q{.} . "$parameters.$identifier.aug";
    $aug_output = safename($aug_output);

    open my $AUG_SCRIPT, '>', $aug_script or die "Can't open AUGUSTUS script $aug_script: $!";

    print $AUG_SCRIPT '#!/bin/bash', "\n\n", ;
    print $AUG_SCRIPT '    augustus --strand=both --genemodel=partial --singlestrand=false --alternatives-from-sampling=true --alternatives-from-evidence=true',
                      ' --minexonintronprob=0.1 --minmeanexonintronprob=0.4 --uniqueGeneId=true --maxtracks=5',
                      ' --protein=on --cds=on --introns=on --start=on --stop=on --codingseq=on',
                      ;

    print $AUG_SCRIPT ' --UTR=on' if $utrs;

    print $AUG_SCRIPT ' --species=', $parameters, ' --outfile=',
                      $aug_output,
                      ' --hintsfile=',
                      $cDNA_hint_file,
                      ' --extrinsicCfgFile=/home/schwarz/src/augustus.2.6.1/config/extrinsic/extrinsic.ME.cfg --progress=true',
                      " $genomic_file ;\n",
                      ; 


    print $AUG_SCRIPT '    program_done_e-ping.pl -p done_', "$aug_script ;\n", ;
    print $AUG_SCRIPT "\n";
    close $AUG_SCRIPT or die "Can't close filehandle to AUGUSTUS script $aug_script: $!"; 

    system "chmod +x $aug_script";

    print $MASTER_SCRIPT '    nohup ./', "$aug_script 1>$aug_script.nohup.out.txt 2>$aug_script.nohup.err.txt &\n"
}

print $MASTER_SCRIPT "\n";
close $MASTER_SCRIPT or die "Can't close filehandle to master script $master_script: $!";

system "chmod +x $master_script";


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


