#!/usr/bin/env perl

use strict;
use warnings;
use strict;

my %strain2species = (
    PB2801 => 'brenneri',
    PS1010 => 'angaria',
    JU1422 => 'Csp9_nigoni',
    DF5081 => 'japonica',
    RS2333 => 'P_pacificus',
    PB4641 => 'remanei',
);

my %strain2output = ();

my $header = '#!/bin/bash' . "\n\n";
my $footer = q{};

my @strains = keys %strain2species;
foreach my $strain (@strains) { 
    if (! $strain2species{$strain} ) { 
        die "Questionable strain-species pair: $strain, $strain2species{$strain}\n";
    }
    my $species = $strain2species{$strain};
    my $output  = "$species.$strain.PacBio_31aug2014.subreads.fa";
    $output     = safename($output);
    $strain2output{$strain} = $output;
}

while (my $input = <>) {
    chomp $input;
    # Sample input:
    # /mnt/home/emsch/work/Ngen/pacbio/meyer/DF5081/E01_1/Analysis_Results/m120622_115256_42134_c100359792550000001523028310041207_s1_p0.subreads.fasta
    if ( $input =~ /\A \/mnt\/home\/emsch\/work\/Ngen\/pacbio\/meyer\/ ([^\s\/]+) \/ [^\s\/]+ \/ Analysis_Results \/ [^\s\/]+ \.subreads\.fasta \z/xms ) { 
        my $strain = $1;
        if (! exists  $strain2output{$strain} ) { 
            die "Cannot identify relevant output file for strain $strain, in input: $input\n";
        }
        print $header if $header;
        $header = q{};
        $footer = "\n";
        print "    cat $input >> $strain2output{$strain};\n";
    }
    else { 
        die "This hard-coded jobmaker cannot parse input: $input;\n";
    }
}

print $footer;


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


