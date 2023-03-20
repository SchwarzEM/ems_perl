#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $seg      = q{};
my $gene2cds = q{};
my $proteome = q{};

# These two variables need to be initialized *outside* the $input loop for $proteome:
my $gene_p = q{};
my $tx_p   = q{};

my $data_ref;
my $verbose;
my $help;

GetOptions ( 'seg=s'      => \$seg,
             'gene2cds=s' => \$gene2cds,
             'proteome=s' => \$proteome,
             'verbose'    => \$verbose,
             'help'       => \$help, );

if ( $help or (! $seg) or (! $proteome) ) { 
    die "Format: summarize_psegs_19mar2023.pl\n",
        "    --seg|-s       [FASTA of low-complexity pseg domains]\n",
        "    --gene2cds|-g  [Gene to CDS TSV table]\n",
        "    --proteome|-p  [proteome; both -s and -p assume AUGUSTUS gene/trx names]\n",
        "    --verbose|-b   [print *all* isoforms values for SEG, rather than default of the highest %]\n",
        "    --help|-h      [print this message]\n",
        ;
}

open my $GENE2CDS, '<', $gene2cds;
while (my $input = <$GENE2CDS>) {
    my $gene = q{};
    my $tx   = q{};
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        $gene = $1;
        $tx   = $2;
    }
    else {
        die "From gene-to-CDS table $gene2cds, can't parse: $input\n";
    }
    $data_ref->{'tx'}->{$tx}->{'gene'} = $gene ;
}
close $GENE2CDS;

open my $PROTEOME, '<', $proteome;
while (my $input = <$PROTEOME>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (\S+) /xms ) {
            $tx_p   = $1;
            $gene_p = q{};
            if ( exists $data_ref->{'tx'}->{$tx_p}->{'gene'} ) {
                $gene_p = $data_ref->{'tx'}->{$tx_p}->{'gene'};
            }
            else {
                die "Cannot map transcript/CDS $tx_p to a gene\n";
            }
        }
        else { 
            die "From proteome $proteome, can't parse: $input\n";
        }
    }
    else {
        $input =~ s/\s//g;
        my $len = length($input);
        $data_ref->{'gene'}->{$gene_p}->{'tx'}->{$tx_p}->{'len'} += $len;
    }
}
close $PROTEOME;

open my $PSEG, '<', $seg;
while (my $input = <$PSEG>) { 
    my $gene = q{};
    my $tx   = q{};
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        # Sample input:
        # >HCON_00000030-00001(758-772) transcript=HCON_00000030-00001 gene=HCON_00000030
        # >HCON_00022890-00001(6-76) transcript=HCON_00022890-00001 gene=HCON_00022890
        if ( $input =~ /\A [>] (\S+) \( (\d+) \- (\d+) \) \s /xms ) { 
            $tx       = $1;
            my $start = $2;
            my $stop  = $3;

            if ( exists $data_ref->{'tx'}->{$tx}->{'gene'} ) {
                $gene = $data_ref->{'tx'}->{$tx}->{'gene'};
            }
            else {
                die "Cannot map transcript/CDS $tx to a gene\n";
            }

            my $seg_span = $stop - $start + 1;
            $data_ref->{'gene'}->{$gene}->{'tx'}->{$tx}->{'seg_span'} += $seg_span;
        }
        else { 
            die "From pseg proteome $seg, can't parse: $input\n";
        }
    }
}
close $PSEG;

my @genes = grep { /\S/ } sort keys %{ $data_ref->{'gene'} };
my $header = "Gene\tPsegs\n";
foreach my $gene1 (@genes) { 
    my $annots_ref;
    my @annots = ();
    my @txs    = grep { /\S/ } sort keys %{ $data_ref->{'gene'}->{$gene1}->{'tx'} };
    foreach my $tx1 (@txs) { 
        if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'}) { 
            my $length   = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'len'};
            my $seg_span = 0;
            if ( exists $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'}) { 
                $seg_span = $data_ref->{'gene'}->{$gene1}->{'tx'}->{$tx1}->{'seg_span'};
            }
            if ( $seg_span > 0 ) { 
                my $orig_fraction = ($seg_span/$length);
                my $fraction = $orig_fraction;
                $fraction = sprintf "%.2f", $fraction;
                my $annot = "$fraction ($seg_span/$length)";
                $annots_ref->{'annot'}->{$annot} = $orig_fraction;
            }
        }
        else { 
            die "Failed to get length for protein product of $tx1\n";
        }
    }
    my $annot_text = q{};
    if (exists $annots_ref->{'annot'}) {
        @annots = grep { /\S/ } sort keys %{ $annots_ref->{'annot'} };
        if ($verbose) {  
            $annot_text = join '; ', @annots;
        }
        else { 
            $annot_text = $annots[-1];        
        }

        # Note that I earlier used:
        #     @annots = sort { $annots_ref->{'annot'}->{$b} <=> $annots_ref->{'annot'}->{$a} } grep { /\S/ } sort keys %{ $annots_ref->{'annot'} };
        #     $annot_text = $annots[0];
        # for the non-verbose option.  But it turns out that a straight ASCII-betical sort of the texts puts my preferred summary at the end anyway!
        # I.e., ASCII-betical picks the largest SEG%, and then the largest AA span if SEG% values tie.  Perfect.
    }
    print $header if $header;
    $header = q{};
    if ( $annot_text =~ /\S/xms ) { 
        $annot_text = 'SEG: ' . $annot_text;
    }
    print "$gene1\t$annot_text\n";
}

