#!/usr/bin/env perl

# get_gene_area_table.pl -- Erich Schwarz <emsch@caltech.edu>, 10/20/2011.
# Purpose: given WormMart output, make table of gene 'regions' from 5' flank to most 3'-ward exon, despite flakiness of WormMart itself on what is "intergenic".

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $input_table = q{};

my $wbgene      = q{};
my $pubgene     = q{};  
my $chr         = q{};
my $strand      = q{};
my $fivePstart  = q{};
my $fivePend    = q{};
my $threePstart = q{};
my $threePend   = q{};
my $seq_ident   = q{};
my $block_start = q{};
my $flank_end   = q{};
my $block_end   = q{};
my $max_flank   = 0;
my $warn_text   = q{};

# Pretty human-readable and normal terms for Watson vs. Crick strands:
my %strand2text = (  1 => '(+)', 
                    -1 => '(-)', );

# If I recall correctly, I previously found that getting Perl to parse '+' or '-' was 
#      way more trouble than it was worth -- hence this peculiar jargon that I put into TSV lines...
my %strand2script_parse = (  1 => 'FOR',
                            -1 => 'REV', );

my $data_ref;
my $help;

GetOptions ( 'table=s' => \$input_table,
             'max=i'   => \$max_flank,
             'help'    => \$help,         );

if ( $help or (! $input_table) ) { 
    die "Format: get_gene_area_table.pl\n",
        "    --table|-t    [WormMart table from which gene regions coordinates are to be extracted]\n",
        "    --max|-m      [Optionally: maximum size of 5' flank in nt, starting from 3'-most 5' flank residue]\n",
        "    --help|-h     [print this message]\n",
        ;
}

if ($max_flank) { 
    if ( (! looks_like_number($max_flank) ) or ( $max_flank < 0 ) or ( $max_flank != int($max_flank) ) ) { 
        die "If a maximum flank size is specified, it should be a nonnegative integer, not $max_flank!\n";
    }
}

# Note that this script was written for a WS220 WormMart data table with a very arbitrary field order.
# Specifically, the data fields were:
#
# Gene WB ID      Gene Public Name        Sequence Name (Transcript)      Transcript Type Chr Name        
#	Strand  Start (bp)      End (bp)        Exon Rank       Exon Start (bp) Exon End (bp)   
#	5' Intergenic Start (bp)        5' Intergenic End (bp)  3' Intergenic Start (bp)        3' Intergenic End (bp)

open my $TABLE, '<', $input_table or die "Can't open input WormMart table $input_table: $!";
while (my $input = <$TABLE>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t      # $wbgene
                       (\S+) \t            # $pubgene
                       \S+ \t

                       # 'coding' is required
                       coding \t 

                       (\S+) \t             # $chr
                       (\S+) \t             # $strand
                       (?: [^\t]* \t){2}

                       # Exon rank 1 is required
                       [1] \t 

                       # For this, exon_start and exon_end are not needed!
                       \d+ \t \d+ \t 

                       # These four points are all we need, and we can even enforce sane order later:
                       (\d+) \t (\d+) \t     # $fivePstart, $fivePend
                       (\d+) \t (\d+)        # $threePstart, $threePend
                    /xms ) { 

        $wbgene      = $1;
        $pubgene     = $2;
        $chr         = $3;
        $strand      = $4;
        $fivePstart  = $5; 
        $fivePend    = $6;
        $threePstart = $7;
        $threePend   = $8; 

        # Make a single name for the block; enforce uniqueness.
        $seq_ident = $wbgene . q{|} . $pubgene;
        if ( exists $data_ref->{'seq_ident'}->{$seq_ident} ) { 
            die "Redundant seq. ident. $seq_ident in: $input\n";
        }

        # Enforce clear choice of strand:
        if ( ( $strand != 1 ) and ( $strand != -1 ) ) { 
            die "Can't parse strandedness of: $input\n";
        }

        # Compute what the coordinates should be, given strand, exon end, and 5' neighbor end.
        # Usually what I'm doing here makes zero difference, because the gene is sensibly annotated.
        # In some rare cases the nt annotations are screwy because of genes-within-genes or such problems.
        # In those case, what follows will *enforce* a commonsensical interpretation!  It may or may not be 
        #    ideal, or even correct, but it will at least give a chance of sanity.

        my @putative_gene_coords = ($fivePstart, $fivePend, $threePstart, $threePend);
        my @gene_boundary_coords = ();

        if ( $strand == 1 ) { 
            @gene_boundary_coords = sort { $a <=> $b } @putative_gene_coords;
        }
        if ( $strand == -1 ) {
            @gene_boundary_coords = reverse sort { $a <=> $b } @putative_gene_coords;
        }

        # Now, get commonsensical coords of gene:
        $block_start = $gene_boundary_coords[0];
        $flank_end   = $gene_boundary_coords[1];
        $block_end   = $gene_boundary_coords[2];

        # Also, for *every* gene we care about at all, record the edges of the gene in the genome,
        #    so that we can subsequently absolutely ensure that dippy WormMart doesn't give us 5' flanks
        #    which actually wander through kilobases of opposite-strand genes -- something it seems to 
        #    want to do, left to its own devices.

        # Use, where appropriate, either the "+ $strand" or the "- $strand" correction.
        #     Partly this is to but the gene points inside an exon terminus;
        #     but mainly this is to avoid having the start of a flank 
        #     *instantly* (and prematurely) shut down the EDGE loop below.

        $data_ref->{'gene_border'}->{'chr'}->{$chr}->{'nt'}->{($flank_end + $strand)} = 1;
        $data_ref->{'gene_border'}->{'chr'}->{$chr}->{'nt'}->{($block_end - $strand)} = 1;

        # If we've specified a maximum 5' flank size, then trim anything larger down to our chosen limit:
        if ($max_flank) { 
            if ( ( $strand == 1 )  and ( $block_start < ( $flank_end - $max_flank + 1 ) ) ) {
                $block_start = ( $flank_end - $max_flank + 1 );
            }
            if ( ( $strand == -1 ) and ( $block_start > ( $flank_end + $max_flank - 1 ) ) ) {
                $block_start = ( $flank_end + $max_flank - 1 );
            }
        }

        # For either value of $strand, this command should put us on the last nt of the last exon:
        $block_end   = $block_end - $strand;

        # Make a header text for the file:
        my $header = $seq_ident . q{   [} . $chr . q{:} . $block_start . q{-} . $block_end . "; $strand2text{$strand} strand]";

        # Store data.  Note that no attempt is made to avoid overwrites where one gene gets 2+ lines,
        #     because the WormMart table *has* such lines and I want robustness here.

        $data_ref->{'seq_id'}->{$seq_ident}->{'chr'}          = $chr;
        $data_ref->{'seq_id'}->{$seq_ident}->{'block_start'}  = $block_start;
        $data_ref->{'seq_id'}->{$seq_ident}->{'flank_end'}    = $flank_end;
        $data_ref->{'seq_id'}->{$seq_ident}->{'block_end'}    = $block_end;
        $data_ref->{'seq_id'}->{$seq_ident}->{'strand'}       = $strand;
    }
}
close $TABLE or die "Can't close filehandle to input WormMart table $input_table: $!";

foreach my $seq_id1 (sort keys %{ $data_ref->{'seq_id'} } ) { 
    # Recover naive values:
    $chr         = $data_ref->{'seq_id'}->{$seq_id1}->{'chr'};
    $block_start = $data_ref->{'seq_id'}->{$seq_id1}->{'block_start'};
    $flank_end   = $data_ref->{'seq_id'}->{$seq_id1}->{'flank_end'};
    $block_end   = $data_ref->{'seq_id'}->{$seq_id1}->{'block_end'};
    $strand      = $data_ref->{'seq_id'}->{$seq_id1}->{'strand'};

    # Then, check if the 'flank' actually overruns genes:
    my @outward_ends   = ($flank_end, $block_start);
    my @outward_coords = ();
    if ( $flank_end > $block_start ) { 
        @outward_ends   = reverse @outward_ends;
        @outward_coords = ($outward_ends[0]..$outward_ends[1]);
        @outward_coords = reverse @outward_coords;
    }
    if ( $flank_end <= $block_start ) { 
        @outward_coords = ($outward_ends[0]..$outward_ends[1]);
    }
    EDGE: foreach my $nt (@outward_coords) { 
        if ( exists $data_ref->{'gene_border'}->{'chr'}->{$chr}->{'nt'}->{$nt} ) { 
            $block_start = $nt + $strand;
            last EDGE;
        }
    }
    # Once I have got the final values for everything, make a header line:
    my $header = $seq_id1 . q{   [} . $chr . q{:} . $block_start . q{-} . $block_end . "; $strand2text{$strand} strand]";

    # Finally, export the data:
    print "$chr\t$block_start\t$block_end\t$header\t$strand2script_parse{$strand}\n";
}

