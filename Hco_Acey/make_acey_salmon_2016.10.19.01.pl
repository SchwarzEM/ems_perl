#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Cwd;
use List::MoreUtils qw(uniq);

my $infile      = q{};
my $curr_dir    = getcwd;
my $i           = 0;
my %file2script = ();

my @prefixes = ();

my $data_ref;

my $help;

GetOptions ( 'infiles=s' => \$infile,
             'help'      => \$help,   );

if ( $help or (! $infile) ) {
    die "Format: make_acey_salmon_2016.10.19.01.pl\n",
        "    --infile|-i   <input table, with prefix-file(s) pairs>\n",
        "    --help|-h     [print this message]\n",
        ;
}

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $prefix     = $1;
        my $input_data = $2;

        if ( $input_data !~ / \.(?:fasta|fq) (?:\.gz){0,1} \z/xms ) {
            die "Cannot accept data as a valid FastQ file: $input_data\n";
        }

        $i++;
        my $j = sprintf "%02u", $i;

        my $output_script = "job_acey_salmon_2016.10.19.$j.sh";
        $output_script    = safename($output_script);

        my $output_prefix = $prefix . ".salmon.2016.10.19.$j";

        if ( exists $data_ref->{'prefix'}->{$prefix} ) {
            die "Redundant input prefix: $prefix\n";
        }

        push @prefixes, $prefix;

        $data_ref->{'prefix'}->{$prefix}->{'input_data'}    = $input_data;
        $data_ref->{'prefix'}->{$prefix}->{'output_script'} = $output_script; 
        $data_ref->{'prefix'}->{$prefix}->{'output_prefix'} = $output_prefix;
    }
}    

$i = 1;

foreach my $prefix (@prefixes) {
    my $input_data    = $data_ref->{'prefix'}->{$prefix}->{'input_data'};
    my $output_script = $data_ref->{'prefix'}->{$prefix}->{'output_script'};
    my $output_prefix = $data_ref->{'prefix'}->{$prefix}->{'output_prefix'};

    open my $OUT, '>', $output_script;

    print $OUT '#!/bin/bash -login', "\n";
    print $OUT '#PBS -l walltime=001:00:00', "\n";
    print $OUT '#PBS -l nodes=1:ppn=8', "\n";
    print $OUT '#PBS -l mem=32gb', "\n";
    print $OUT "#PBS -N $output_script\n";
    print $OUT '#PBS -q main', "\n";
    print $OUT '#PBS -M ems394@cornell.edu', "\n";
    print $OUT '#PBS -m abe', "\n";
    print $OUT '#PBS -A ged', "\n";
    print $OUT '#PBS -r n', "\n";
    print $OUT '#PBS -V', "\n";
    print $OUT "cd $curr_dir ;\n";
    print $OUT 'export PATH="$PATH:/mnt/home/emsch/src/Salmon-0.7.2_linux_x86_64/bin" ;', "\n";

    print $OUT "salmon --no-version-check quant --threads 8 ";
    print $OUT "--index  /mnt/home/emsch/work/Acey/immuno_rnaseq/salmon/dbs/ancylostoma_ceylanicum.PRJNA231479.WBPS6.CDS_transcripts ";
    print $OUT "--libType A ";
    print $OUT "--unmatedReads $input_data ";
    print $OUT "--output $output_prefix ";
    print $OUT "--seqBias --gcBias ";
    print $OUT "--geneMap /mnt/home/emsch/work/Acey/immuno_rnaseq/salmon/dbs/Acey_PRJNA231479_WBPS6.cds2gene.tsv.txt ";
    print $OUT "--numBootstraps 100 ;";

    print $OUT "\n";

    if ( exists $prefixes[$i] ) {
        my $next_output_script = $data_ref->{'prefix'}->{$prefixes[$i]}->{'output_script'};
        print $OUT "qsub $next_output_script ;\n";
    }

    close $OUT;
    $i++;
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

