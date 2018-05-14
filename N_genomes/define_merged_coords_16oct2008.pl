#!/usr/bin/env perl

# define_merged_coords_16oct2008.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/16/2008.  LEGACY version, kept for ability to reproduce older work.
# Purpose: given a table of coords in elegans genome, define merged coords.
# Where possible, preserve strand orientations (may not always be able to do this).

use strict;
use warnings;

my $marked_nt_ref;
my $block_list_ref;

my %chrs = ();
my @chr_names = qw( I II III IV V X MtDNA );
foreach my $chr_name (@chr_names) { 
    $chrs{$chr_name} = 1;
}

while (my $input = <>) { 
    my $chr  = q{};
    my $nt1  = q{};
    my $nt2  = q{};
    my $ori  = q{};
    my $gene = q{};
    chomp $input;
    if ( $input =~ / \A (\S+) \t (\d+) \t (\d+) \t (\S+) \t (\S+)/xms ) { 
        $chr  = $1;
        $nt1  = $2;
        $nt2  = $3;
        $ori  = $4;
        $gene = $5;
        if (! exists $chrs{$chr} ) { 
            warn "Can't parse chromosome: $input\n";
        }
        if ( $nt1 > $nt2) { 
            ($nt1, $nt2) = ($nt2, $nt1);
        }
        foreach my $i ($nt1..$nt2) { 
            $marked_nt_ref->{$chr}->{$i}->{'gene'}->{$gene} = 1;
            $marked_nt_ref->{$chr}->{$i}->{'ori'}->{$ori} = 1;
        }
    }
}

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
              $gene_text,
              "\t",
              $orient_text,
              "\n",
              ;
    }
}

