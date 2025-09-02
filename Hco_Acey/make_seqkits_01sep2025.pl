#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list = q{};
$list    = $ARGV[0] if $ARGV[0];

if (! $list ) {
    die "Format: make_seqkits_01sep2025.pl [list] > [script]\n"
}

open my $LIST, '<', $list;

while (my $input = <$LIST>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\d+) \s+ (\d+) \z/xms ) {
         my $stem    = $1;
         my $start   = $2;
         my $stop    = $3;
         my $range_c = "$start:$stop";
         my $range_d = "$start-$stop";
         my $infile  = "$stem.pep.fa";
         my $outfile = "$stem.$range_d.pep.fa";
         print "extract_fasta_subset.pl -f seqs/possible_new_ligands.pep.fa -l $stem > $infile ;\n";
         print "seqkit subseq -R -r $range_c $infile > $outfile ;\n";
    }
    else {
        die "From list $list, cannot parse: $input\n";
    }
}
