#!/usr/bin/env perl

# get_real_flanks.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/15/2008.
# Purpose: given a WormMart table, get non-trivial coordinates of 5' flanks of real 5' UTRs.
# N.B. The basic idea is to get things most likely to be promoters; hence the requirement for a true 5' UTR.

use strict;
use warnings;
use Getopt::Long;

my $length  = q{};
my $ok_type = q{};

GetOptions ( "length=i" => \$length,
             "type=s"   => \$ok_type,   );

# Defaults:
$length  |= 100;
$ok_type |= 'coding';

# WBGene00000003  aat-2   F07C3.7 9246327 9246334 9246327 9248004 V       9244386 9246334 F07C3.7 coding

my $orientation;
my $wb_id;
my $pubname;
my $seqname;
my $utr5start;
my $utr5end;
my $inter5start;
my $inter5end;
my $chr;
my $gene_start;
my $gene_stop;
my $tx_type;

my $segments_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A (WBGene\d+)           # Gene WB ID
                        \t 
                        (\S+)                 # Gene Public Name
                        \t 
                        (\S+)                 # Sequence Name (Gene)
                        \t 
                        (\d+) \t (\d+)        # 5'UTR Start (bp); 5'UTR end (bp)
                        \t 
                        (\d+) \t (\d+)        # 5' Intergenic Start (bp); 5' Intergenic End (bp)
                        \t 
                        (\S+)                 # Chr Name
                        \t (\d+) \t (\d+) \t    # Start (bp); End (bp)
                        \S+                   # Sequence Name (Transcript)
                        \t 
                        (\S+)                 # Transcript Type
                  /xms ) { 
        ( $wb_id,     $pubname, $seqname,
          $utr5start, $utr5end, $inter5start,
          $inter5end, $chr,     $gene_start,
          $gene_stop, $tx_type, ) =  ( $1, $2, $3, 
                                       $4, $5, $6, 
                                       $7, $8, $9, 
                                       $10, $11, );
        $orientation = q{};
        if ( abs($utr5start - $utr5end) > 1 ) { 

            # Note: these criteria are, arguably, naively overconservative.
            # On a set of 9,263 5'-flank sets from WS190, they excluded 1,042.
            # However, they do have the advantage of providing a highly 
            # unambiguous starting set of sequences.

            if ( ( $inter5start < $gene_start ) 
                 and ( $utr5start == $gene_start ) ) { 
                $orientation = 'FOR';
            }
            if ( ( $gene_stop < $inter5end ) 
                 and ( $utr5start > $gene_start ) ) { 
                $orientation = 'REV'; 
            }
            # In WS190, only 49 out of 8,221 transcripts were non-coding:
            if ( ( $orientation ) and ( $tx_type eq $ok_type ) ) { 
                my $select1;
                my $select2;
                if ( $orientation eq 'FOR') { 
                    $select2 = $utr5start - 1;
                    $select1 = $select2 - $length + 1;
                }
                if ( $orientation eq 'REV') {
                    $select1 = $utr5end + 1; 
                    $select2 = $select1 + $length - 1;
                }
                my $gene_id = $wb_id . '|' . $seqname;
                if ( $pubname ne $seqname ) { 
                    $gene_id .= "|$pubname";
                }
                my $segment_coords = "$chr\t$select1\t$select2\t$orientation";
                $segments_ref->{$segment_coords}->{$gene_id} = 1;
            }
        }
    }
}

foreach my $segment (sort keys %{ $segments_ref } ) { 
    foreach my $gene (sort keys %{ $segments_ref->{$segment} } ) { 
        print "$segment\t$gene\n";
    }
}

