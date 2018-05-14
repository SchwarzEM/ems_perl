#!/usr/bin/env perl

# auto_HMerge_genome_split.pl -- Erich Schwarz <ems394@cornell.edu>, 7/14/2014.
# Purpose: automate the job of getting a single assembly split into a larger 10% and a smaller 90% (for HaploMerger matrix-building).

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my @infiles  = ();
my $fraction = 0.1;
my @script   = ();
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'fraction=f'   => \$fraction,
             'help'         => \$help,     );

if ( $help or (! @infiles) ) { 
    die "Format: auto_HMerge_genome_split.pl\n",
        "    --infile|-i    <input stream/files>\n",
        "    --fraction|-f  [decimal fraction, greater than 0.00 and less than 1.00; default is 0.10]\n",
        "    --help|-h      [print this message]\n",
        ;
}

if ( (! looks_like_number($fraction) ) and ( $fraction <= 0 ) or ( $fraction >= 1 ) ) { 
    die "Fraction must be decimal fraction >0, <1\n";
}

foreach my $infile (@infiles) {
    my $basename = basename($infile);
    $basename =~ s/\.fa\z//;

    my $top_file    = $basename . '_top.fa';
    my $bottom_file = $basename . '_bottom.fa';
    my $top_list    = $basename . '.top.list.txt';

    $top_file       = safename($top_file);
    $bottom_file    = safename($bottom_file);
    $top_list       = safename($top_list);

    my $count_nt = `count_fasta_residues.pl -i $infile | grep "Total nt:"`;
    chomp $count_nt;
    if ( $count_nt !~ /\A Total [ ] nt: [ ]+ \S*\d \z/xms ) { 
        die "Didn't manage to parse the sequence size: $count_nt\n";
    }
    if ( $count_nt =~ /\A Total [ ] nt: [ ]+ (\S*\d) \z/xms ) {
        $count_nt = $1;
    }
    $count_nt =~ s/,//g;
    $count_nt = ($count_nt * $fraction);
    $count_nt = int $count_nt;
    if ( $count_nt < 1 ) { 
        die "Count nt ($count_nt) is too small\n";
    }

    my $command1 = "    fasta_size_subset.pl -i $infile -f --max $count_nt > $top_file ;\n";
    push @script, $command1;
    my $command2 = q{    grep '>' } . $top_file . q{ | perl -ne ' m/>(\S+)/; print "$1\n"; ' | sort | uniq > } . " $top_list ;\n";
    push @script, $command2;
    my $command3 = "    extract_fasta_subset.pl -e -l $top_list -f $infile > $bottom_file ;\n";
    push @script, $command3;
    my $command4 = "    nohup gzip -9 $infile &\n";
    push @script, $command4;
    my $command5 = "    nohup gzip -9 $top_file &\n";
    push @script, $command5;
    my $command6 = "    nohup gzip -9 $bottom_file &\n";
    push @script, $command6;
    my $command7 = "\n";
    push @script, $command7;
}

my $header = '#!/bin/bash' . "\n\n";

if (@script) {
    print $header;
    print @script;
}

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
 
