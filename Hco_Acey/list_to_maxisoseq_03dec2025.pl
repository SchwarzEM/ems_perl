#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $list = q{};
$list    = $ARGV[0] if $ARGV[0];

if (! $list ) {
    die "Format: list_to_maxisoseq_03dec2025.pl [seqlist] > [extraction commands]\n";
}

my $bigfasta = '/ocean/projects/mcb190015p/schwarze/heligmo/orig_proteomes/nxHelBake1.1.primary.final_annotations.proteins.rev1_wHES.max_isos.fa';

open my $LIST, '<', $list;
while (my $seq = <$LIST>) {
    chomp $seq ;
    my $command = "extract_fasta_subset.pl -w -f $bigfasta -l $seq > $seq.fa ;";
    print "$command\n";
}
close $LIST;

