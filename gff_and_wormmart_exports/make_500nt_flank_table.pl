#!/usr/bin/env perl

# make_500nt_flank_table.pl -- Erich Schwarz <emsch@caltech.edu>, 10/16/2011.
# Purpose: given WormMart output, make nonredundant table of 500-nt flanks, constrained by boundaries of adjacent gene.  Make no effort to prune to boundary of adjacent gene (this works much better if exons are masked than if they're not!).

use strict;
use warnings;

my $wbgene      = q{};
my $pubgene     = q{};  
my $tx          = q{};   
my $chr         = q{};
my $strand      = q{};
my $exon_start  = q{};
my $exon_end    = q{};
my $seq_ident   = q{};
my $block_start = q{};
my $block_end   = q{};
my $warn_text   = q{};

# Pretty human-readable and normal terms for Watson vs. Crick strands:
my %strand2text = (  1 => '(+)', 
                    -1 => '(-)', );

# If I recall correctly, I previously found that getting Perl to parse '+' or '-' was 
#      way more trouble than it was worth -- hence this peculiar jargon that I put into TSV lines...
my %strand2script_read = (  1 => 'FOR',
                           -1 => 'REV', );

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t 
                       (\S+) \t 
                       (\S+) \t

                       # 'coding' is required
                       coding \t 

                       (\S+) \t
                       (\S+) \t
                       (?: [^\t]* \t){2}

                       # Exon rank 1 is required
                       [1] \t 

                       # $exon_start, $exon_end
                       (\d+) \t (\d+) \t 
                    /xms ) { 

        $wbgene      = $1;
        $pubgene     = $2;
        $tx          = $3;
        $chr         = $4;
        $strand      = $5;
        $exon_start  = $6;
        $exon_end    = $7;

        # Make a single name for the block; enforce uniqueness.
        $seq_ident = $wbgene . q{|} . $pubgene . q{|}. $tx;
        if ( exists $data_ref->{'seq_ident'}->{$seq_ident} ) { 
            die "Redundant seq. ident. $seq_ident in: $input\n";
        }

        # Compute what the coordinates should be, given strand, exon end, and 5' neighbor end:
        if ( $strand == 1 ) { 
            $block_end   = $exon_start - 1;
            $block_start = $exon_start - 500;
        }
        if ( $strand == -1 ) {
            $block_start = ($exon_end + 1);
            $block_end   = ($exon_end + 500);
        }

        # Enforce clear choice of strand:
        if ( ( $strand != 1 ) and ( $strand != -1 ) ) { 
            die "Can't parse strandedness of: $input\n";
        }

        # Make a header text for the file:
        my $header = $seq_ident . q{   [} . $chr . q{:} . $block_start . q{-} . $block_end . "; $strand2text{$strand} strand]";

        # Basic sanity checks:
        if ( $block_end < $block_start ) {
            $block_end = $block_start;
            $warn_text = "WARNING: Original block end was lower than block start $block_start for: $input";
            warn "$warn_text\n";
            $header = $header . "  [$warn_text]";
        }
        if ( $block_start < 1 ) {
            $block_start = 1;
            $warn_text = "WARNING: Original block start was lower than 1 for: $input";
            warn "$warn_text\n";
            $header = $header . "  [$warn_text]";
        }
        if ( $block_end < 1 ) {
            $block_end = 1;
            $warn_text = "WARNING: Original block end was lower than 1 for: $input";
            warn "$warn_text\n";
            $header = $header . "  ($warn_text)";
        }

        # Export data:
        print "$chr\t$block_start\t$block_end\t$header\t$strand2script_read{$strand}\n";
    }
}

