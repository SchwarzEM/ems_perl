#!/usr/bin/env perl

# merge_nonelegans_PFAM_lists.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/6/2010.
# Purpose: for working up PS1010 paper table -- merge several partial tables into one.

use strict;
use warnings;

my $motif_txt = q{};
my $species   = q{};

my $data_ref;
my @motif_data = ();
my @motif_spp  = ();

while (my $input = <>) { 
    chomp $input;

    # Sample input line:
    # 
    # 6	PF01842.18	ACT	ACT domain	"P. pacificus"

    if ( $input =~ / \A ( \d+ \t [^\t]+ \t [^\t]+ \t [^\t]+ ) \t \" ( [^\"\t]+ ) \" /xms ) { 
        $motif_txt = $1;
        $species  = $2;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
    $data_ref->{$motif_txt}->{$species} = 1;
}

# Impose a gene-count sort on *top* of an ASCIIbetal plain 'sort'.
@motif_data = sort { &get_motif_count($b) <=> &get_motif_count($a) } sort keys %{ $data_ref };

foreach my $motif_datum (@motif_data) { 
    @motif_spp = sort keys %{ $data_ref->{$motif_datum} };
    my $spp_list = join ", ", @motif_spp;
    print "$motif_datum\t$spp_list\n";
}

sub get_motif_count { 
    my $_input = $_[0];
    my $_count = 0;
    if ( $_input =~ /\A (\d+) \t /xms ) { 
        $_count = $1;
    }
    else { 
       die "Can't parse motif count in: $_input\n";
    }
    return $_count;
}

