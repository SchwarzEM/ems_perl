#!/usr/bin/env perl

# give_wormpep_full_gene_names.pl -- Erich Schwarz <ems394@cornell.edu>, 6/24/2018
# Purpose: given wormpep with 1-protein-per-gene, rename the proteins so that they have full gene names (useful for getting BlastP reports that are fully informative about the gene being hit).

use strict;
use warnings;
use autodie;

my $wormpep = q{};
$wormpep    = $ARGV[0] if $ARGV[0];

my %cds2fullname = ();
my %possible_cds = ();
my %seen         = ();

# Typical wormpep190 lines:
# >4R79.1b        CE39659 WBGene00003525  locus:nas-6 ..
# >AC7.3  CE07653 WBGene00014997

if (! -e $wormpep) {
    die "Format: give_wormpep_full_gene_names.pl [wormpep subset, one-protein-per-gene] > [full-gene-renamed wormpep]\n";
}

open my $WORMPEP, '<', $wormpep;
while (my $input = <$WORMPEP>) { 
    chomp $input;

    if ($input =~ /\A [>] ( (\S+) \b .+ (WBGene\d+) .*) \z/xms) { 
        my $header      = $1;
        my $cds         = $2;
        my $gene        = $3;
        my $locus       = q{};

        $cds =~ s/[a-z]\z//;

        my @id_tags = ($gene, $cds);

        if ( $input =~ / locus = (\S+) /xms ) {
            $locus = $1;
            push @id_tags, $locus;
        }

        my $full_id = join '|', @id_tags;

        if ( exists $seen{$full_id} ) {
            die "Redundant gene name ($full_id) for protein: $input\n";
        }
        else {
            print ">$full_id  $header\n";
            $seen{$full_id} = 1;
        }
    }

    elsif ( $input =~ /\A [>] /xms) {
        die "Cannot parse FASTA header: $input\n";
    }

    else {
        print "$input\n";
    }
}
close $WORMPEP;

