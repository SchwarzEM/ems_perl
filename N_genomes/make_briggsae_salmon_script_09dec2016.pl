#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $geneset_table = q{};
my $read_table    = q{};

my @gene_sets   = ();
my %genes2index = ();
my %genes2txmap = ();
my @commands    = ();
my %seen        = ();

$geneset_table = $ARGV[0] if $ARGV[0];
$read_table    = $ARGV[1] if $ARGV[1];

if ( (! $geneset_table) or (! $read_table) ) {
    die "Format: make_briggsae_salmon_script_09dec2016.pl\n",
        "            [gene set tab-delimited table: gene set, salmon index, tx2gene map]\n",
        "            [read file tab-delimited table: read file, annotation]\n",
        "            > [qsub-able script]\n",
        ;
}

open my $GENESET, '<', $geneset_table;
while (my $input = <$GENESET>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
        my $gene_set     = $1;
        my $salmon_index = $2;
        my $tx2gene_map  = $3;

        if ( ( exists $genes2index{$gene_set} ) or ( exists $genes2txmap{$gene_set} ) ) {
            die "In gene set table $geneset_table, redundant gene set $gene_set (in: $input)\n";
        }

        push @gene_sets, $gene_set;

        $genes2index{$gene_set} = $salmon_index;
        $genes2txmap{$gene_set} = $tx2gene_map;
    }
    else { 
        die "In gene set table $geneset_table, cannot parse text line: $input\n";
    }
}
close $GENESET;

open my $READ, '<', $read_table;
while (my $input = <$READ>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $read_file  = $1;
        my $read_annot = $2;
        my $srr_acc    = q{};

        if (exists $seen{$read_file} ) {
            die "In read file $read_file, redundant read file $read_file (in: $input)\n";
        }
        $seen{$read_file} = 1;

        if (! -r $read_file) {
            die "Cannot read putative readfile: $read_file\n";
        }

        if ( $read_file =~ /(SRR\d+)/xms ) { 
            $srr_acc = $1;
        }
        else { 
            die "Cannot identify SRR accession in: $read_file\n";            
        }

        foreach my $gene_set (@gene_sets) {
            my $command = "salmon --no-version-check quant --threads 8 --index $genes2index{$gene_set} "
                          . "--libType A --unmatedReads $read_file --output $gene_set.$srr_acc.$read_annot "
                          . "--seqBias --gcBias --geneMap $genes2txmap{$gene_set} --numBootstraps 100 ;"
                          ;
            push @commands, $command;
        }
    }
    else {
        die "In read table $read_table, cannot parse text: $input\n";
    }
}
close $READ;

@commands  = uniq(@commands);

print '#!/bin/bash -login', "\n";
print '#PBS -l walltime=024:00:00', "\n";
print '#PBS -l nodes=1:ppn=8', "\n";
print '#PBS -l mem=32gb', "\n";
print '#PBS -N job_CURRENT.sh', "\n";
print '#PBS -q main', "\n";
print '#PBS -M ems394@cornell.edu', "\n";
print '#PBS -m abe', "\n";
print '#PBS -A ged', "\n";
print '#PBS -r n', "\n";
print '#PBS -V', "\n";
print 'cd WORKING_DIRECTORY ;', "\n";
print 'export PATH="$PATH:/mnt/home/emsch/src/SalmonBeta-0.7.0_linux_x86_64/bin" ;', "\n";

foreach my $command (@commands) {
    print "$command\n";
}


