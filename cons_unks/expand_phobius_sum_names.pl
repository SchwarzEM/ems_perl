#!/usr/bin/env perl

use strict;
use warnings;
use strict;

my $fasta = q{};
my $summ  = q{};

if ( (! $ARGV[0]) or (! $ARGV[1]) ) {
    die "Format: expand_phobius_sum_names.pl [FlyPep, 2016 version] [FlyPep-based Phobius summary]\n";
}
else {
    $fasta = $ARGV[0];
    $summ  = $ARGV[1];
}

my $data_ref;

open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ /\A [>] \S+ \s .* name [=] (\S+) [-] P [A-Z]+ [;] \s+ parent [=] (FBgn\d+) [,] /xms ) { 
        my $name    = $1;
        my $gene_id = $2;
        if ( ( exists $data_ref->{'gene_id'}->{$gene_id} ) and ( $data_ref->{'gene_id'}->{$gene_id} ne $name ) ) {
            die "Contradictory names for gene ID $gene_id: \"$data_ref->{'gene_id'}->{$gene_id}\" vs. \"$name\"\n";
        }
        $data_ref->{'gene_id'}->{$gene_id} = $name ;
    }
    elsif ( $input =~ /\A [>] \S /xms ) { 
        die "Cannot parse FASTA header line: $input\n";
    }
}
close $FASTA;

open my $SUMM, '<', $summ;
while (my $input = <$SUMM>) {
    chomp $input;
    my $output = $input;
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) { 
        my $gene_id = $1;
        my $annot   = $2;
        if ( exists $data_ref->{'gene_id'}->{$gene_id} ) {
            $gene_id = $gene_id . '|' . $data_ref->{'gene_id'}->{$gene_id} ;
            $output = "$gene_id\t$annot";
        }
        elsif ( $gene_id ne 'Gene' ) { 
            die "Cannot convert gene name in Phobius summary input line: $input\n";
        }
    }
    print "$output\n";
}

