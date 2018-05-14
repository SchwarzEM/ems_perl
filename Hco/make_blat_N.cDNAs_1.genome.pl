#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my $master_script = q{};
my @cDNA_files    = ();
my $genomic_file  = q{};
my $blat_script   = q{};

my $help;

GetOptions ( 'master=s'  => \$master_script,
             'cdna=s{,}' => \@cDNA_files,
             'genome=s'  => \$genomic_file,
             'help'      => \$help, );

if ( $help or (! $master_script) or (! @cDNA_files) or (! $genomic_file) ) {
    die "\n",
        "Format: make_blat_N.cDNAs_1.genome.pl\n",
        "        --master|-m  [name of intended master script file]\n",
        "        --cdna|-a    [name(s) of one or more existing cDNA files]\n",
        "        --genome|-g  [name of one existing genomic DNA file]\n",
        "        --help|-h    [print this message]\n",
        "\n",
        ;
}

$master_script = safename($master_script);

open my $MASTER_SCRIPT, '>', $master_script or die "Can't open master script $master_script: $!";

print $MASTER_SCRIPT '#!/bin/bash', "\n\n", ;

foreach my $cDNA (@cDNA_files) { 
    my $cDNA_base = basename $cDNA;
    my $gen_base  = basename $genomic_file;

    $blat_script = 'job_blat_' . $cDNA_base . '_vs_' . $gen_base . '.sh';
    $blat_script = safename($blat_script);

    my $psl_file = $gen_base . '_vs_' . $cDNA . '.psl';
    $psl_file = safename($psl_file);

    open my $BLAT_SCRIPT, '>', $blat_script or die "Can't open BLAT script $blat_script: $!";

    print $BLAT_SCRIPT '#!/bin/bash', "\n\n", ;
    print $BLAT_SCRIPT "    blat -minIdentity=92 $genomic_file $cDNA $psl_file ;\n";
    print $BLAT_SCRIPT '    program_done_e-ping.pl -p done_', "$blat_script ;\n", ;
    print $BLAT_SCRIPT "\n";

    close $BLAT_SCRIPT or die "Can't close filehandle to BLAT script $blat_script: $!"; 

    system "chmod +x $blat_script";

    print $MASTER_SCRIPT '    nohup ./', "$blat_script 1>$blat_script.nohup.out.txt 2>$blat_script.nohup.err.txt &\n"
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


