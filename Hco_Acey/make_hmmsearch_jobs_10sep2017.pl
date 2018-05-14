#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Basename;

my @families  = ();
my @proteomes = ();
my $header    = '#!/bin/bash';
my $help;

GetOptions ( 'families=s{,}'  => \@families,
             'proteomes=s{,}' => \@proteomes,
             'help'           => \$help,   );

if ( $help or (! @families) or (! @proteomes) ) { 
    die "Format: make_hmmsearch_jobs_10sep2017.pl\n",
        "    --families|-f   [HMMs of protein families to use with hmmsearch\n",
        "    --proteomes|-p  [proteomes to use with with hmmsearch\n",
        "    --help|-h       [print this message]\n",
        ;
}

foreach my $proteome (@proteomes) {
    print "$header\n\n" if $header;
    $header = q{};

    my $prot_base = basename($proteome);
    $prot_base =~ s/\.fasta\z//;
    $prot_base =~ s/\.fa\z//;

    foreach my $hmm (@families) {
        my $hmm_base = basename($hmm);
        $hmm_base =~ s/\.hmm\z//;
        print "hmmsearch --cpu 8 --domE 0.5 -o $hmm_base.hmmsearch.$prot_base.txt $hmm $proteome ;\n";
    }

    print "\n";
}
