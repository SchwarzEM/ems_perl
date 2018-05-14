#!/usr/bin/env perl

# define_merged_coords.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/16/2008; significantly modified on 11/6/2013.
# Purpose: given a table of coords in a set of DNA sequences (with chromosomes, scaffolds, etc.), define merged coords.
# Where possible, preserve strand orientations (may not always be able to do this).

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $marked_nt_ref;
my $block_list_ref;

my @chr_names = ();

while (my $input = <>) { 
    my $chr  = q{};
    my $nt1  = q{};
    my $nt2  = q{};
    my $ori  = q{};
    my $gene = q{};
    chomp $input;
    if ( $input =~ / \A (\S+) \t (\d+) \t (\d+) \t (\S+) \t (\S+) /xms ) { 
        $chr  = $1;
        $nt1  = $2;
        $nt2  = $3;
        $ori  = $4;
        $gene = $5;
        if ( $nt1 > $nt2) { 
            ($nt1, $nt2) = ($nt2, $nt1);
        }
        foreach my $i ($nt1..$nt2) { 
            $marked_nt_ref->{$chr}->{$i}->{'gene'}->{$gene} = 1;
            $marked_nt_ref->{$chr}->{$i}->{'ori'}->{$ori} = 1;
        }
        push @chr_names, $chr;
    }
}

# Do it this way to keep the order of chr names present in the input file.
@chr_names = uniq @chr_names;

foreach my $chromosome (@chr_names) { 
    my @coords = sort { $a <=> $b } 
                 keys %{ $marked_nt_ref->{$chromosome} };
    my $block_id = q{};
    foreach my $coord (@coords) { 
        my $prev_nt = $coord - 1;
        my $next_nt = $coord + 1;
        if (! exists $marked_nt_ref->{$chromosome}->{$prev_nt} ) {
            $block_id = $chromosome . ':' . $coord;
            $block_list_ref->{$chromosome}->{$block_id}->{'start'} = "$coord";
        }
        if (! exists $marked_nt_ref->{$chromosome}->{$next_nt} ) { 
            $block_list_ref->{$chromosome}->{$block_id}->{'stop'} = "$coord"; 
        }
        my @genes
            = sort keys %{ $marked_nt_ref->{$chromosome}->{$coord}->{'gene'} };
        my @orientations
            = sort keys %{ $marked_nt_ref->{$chromosome}->{$coord}->{'ori'} };
        foreach my $gene (@genes) {
            $block_list_ref->{$chromosome}->{$block_id}->{'gene'}->{$gene} = 1;
        }
        foreach my $ori (@orientations) {
            $block_list_ref->{$chromosome}->{$block_id}->{'ori'}->{$ori} = 1;
        }
    }
}

foreach my $chromosome (@chr_names) { 
    foreach my $block_id (sort keys %{ $block_list_ref->{$chromosome} } ) { 
        my $gene_text   = join '; ',
            sort keys %{ $block_list_ref->{$chromosome}->{$block_id}->{'gene'} };
        my $orient_text = join '; ', 
            sort keys %{ $block_list_ref->{$chromosome}->{$block_id}->{'ori'} };
        print $chromosome,
              "\t",
              $block_list_ref->{$chromosome}->{$block_id}->{'start'},
              "\t",
              $block_list_ref->{$chromosome}->{$block_id}->{'stop'},
              "\t",
              $block_id,
              "\t",
              $orient_text,
              "\t",
              $gene_text,
              "\n",
              ;
    }
}

