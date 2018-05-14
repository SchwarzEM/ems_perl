#!/usr/bin/env perl

# summarize_psegs.pl -- Erich Schwarz <ems394@cornell.edu>, 9/28/2014.
# Purpose: given the segments of proteome identified by pseg, the proteome's pre-masked version, and a gene-to-CDS table, generate a gene-oriented SEG annotation table.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $gene2cds = q{};
my $seg      = q{};
my $proteome = q{};

my $gene = q{};
my $cds  = q{};

my $data_ref;
my $verbose;
my $help;

GetOptions ( 'gene2cds=s' => \$gene2cds,
             'seg=s'      => \$seg,
             'proteome=s' => \$proteome,
             'verbose'    => \$verbose,
             'help'       => \$help, );

if ( $help or (! $gene2cds) or (! $seg) or (! $proteome) ) { 
    die "Format: summarize_psegs.pl\n",
        "    --gene2cds|-g  [Gene-to-CDS table, used to map CDS/protein names in PSEG-masked proteome]\n",
        "    --proteome|-p  [unmasked source proteome]\n",
        "    --seg|-s       [FASTA of low-complexity pseg domains]\n",
        "    --verbose|-b   [print *all* isoforms values for PSEG, rather than default of the highest %]\n",
        "    --help|-h      [print this message]\n",
        ;
}

open my $GENE, '<', $gene2cds;
while (my $input = <$GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        $gene = $1;
        $cds  = $2;
        if ( ( exists $data_ref->{'cds'}->{$cds}->{'gene'} ) and ( $data_ref->{'cds'}->{$cds}->{'gene'} ne $gene ) ) {
            die "CDS $cds is inconsistently mapped to both $data_ref->{'cds'}->{$cds}->{'gene'} and $gene\n";
        }
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else {
        die "From gene-to-CDS name table $gene2cds, can't parse input line: $input\n";
    }
}
close $GENE;

open my $PROTEOME, '<', $proteome;
while (my $input = <$PROTEOME>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        if ( $input =~ /\A > (\S+) /xms ) {
            $cds  = $1;
            $gene = q{};
            if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
                die "In proteome $proteome, cannot map CDS $cds to gene\n";
            }
            $gene = $data_ref->{'cds'}->{$cds}->{'gene'};
        }
        else {
            die "From proteome $proteome, can't parse: $input\n";
        }
    }
    else {
        $input =~ s/\s//g;
        my $len = length($input);
        $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'len'} += $len;
    }
}
close $PROTEOME;

open my $PSEG, '<', $seg;
while (my $input = <$PSEG>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        my $start = q{};
        my $stop  = q{};

        if ( $input =~ /\A > (\S+) \b .* \( (\d+) \- (\d+) \) /xms ) { 
            $cds   = $1;
            $start = $2;
            $stop  = $3;
        }
        else {
            die "From PSEG file $seg, can't parse: $input\n";
        }

        if ( $input =~ / \( (\d+) \- (\d+) \) .* \( (\d+) \- (\d+) \) /xms ) { 
            die "In PSEG file $seg, anomalous header line whose PSEG-marked coordinates cannot be reliably identified: $input\n";
        }

        if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
            die "In PSEG file $seg, cannot map CDS $cds to gene\n";
        }
        $gene = $data_ref->{'cds'}->{$cds}->{'gene'};

        my $seg_span = $stop - $start + 1;
        $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'seg_span'} += $seg_span;
    }
}
close $PSEG;

my @genes = grep { /\S/ } sort keys %{ $data_ref->{'gene'} };
my $header = "Gene\tPsegs\n";
foreach my $gene1 (@genes) { 
    my $annots_ref;
    my @annots = ();
    my @cdses  = grep { /\S/ } sort keys %{ $data_ref->{'gene'}->{$gene1}->{'cds'} };
    foreach my $cds1 (@cdses) { 
        if ( exists $data_ref->{'gene'}->{$gene1}->{'cds'}->{$cds1}->{'len'}) { 
            my $length   = $data_ref->{'gene'}->{$gene1}->{'cds'}->{$cds1}->{'len'};
            my $seg_span = 0;
            if ( exists $data_ref->{'gene'}->{$gene1}->{'cds'}->{$cds1}->{'seg_span'}) { 
                $seg_span = $data_ref->{'gene'}->{$gene1}->{'cds'}->{$cds1}->{'seg_span'};
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
            die "Failed to get length for protein product of $cds1\n";
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

