#!/usr/bin/env perl

# fasta_length.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/29/2008.
# Purpose: censors too short (e.g., zero-length) or too long FASTAs.

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $prot_name = q{};
my $gene_name = q{};
my $gene2prot_ref;
my %prot2gene = ();
my %prot2len = ();

# Record the incoming FASTA data.

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A > (\S+) .* (WBGene\d+) /xms) { 
        $prot_name = $1;
        $gene_name = $2; 
        $prot2gene{$prot_name} = $gene_name;
        $prot2len{$prot_name} = 0;
    }
    elsif ( $input =~ /\S/xms ) { 
        $prot2len{$prot_name} += length($input);
    }
}

foreach my $prot_name2 (sort keys %prot2gene) { 
    my $gene_name2 = $prot2gene{$prot_name2};
    push @{ $gene2prot_ref->{$gene_name2}->{'lengths'} }, $prot2len{$prot_name2};
}

foreach my $gene_name3 ( sort keys %{ $gene2prot_ref } ) { 
    my @prot_lens = @{ $gene2prot_ref->{$gene_name3}->{'lengths'} };
    my $min_protlen = min @prot_lens;
    my $max_protlen = max @prot_lens;
    print "$gene_name3\t";
    if ( $min_protlen < $max_protlen ) { 
        print "\"$min_protlen-$max_protlen aa\"\n";
    }
    else { 
        print "\"$max_protlen aa\"\n";
    }
}

