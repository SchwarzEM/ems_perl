#!/usr/bin/env perl

# cds2gene_annot.pl -- Erich Schwarz <ems394@cornell.edu>, 3/4/2023.
# Purpose: given a gene2cds table and 1+ columns of annots, make a table of genes and annotations (e.g., GO terms) joined with '; '.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $gene2cds = q{};
my $annot    = q{};
my $title    = q{};

my $no_split;
my $help;

GetOptions ( 'gene2cds=s' => \$gene2cds,
             'annot=s'    => \$annot,
             'title=s'    => \$title,
             'no_split'   => \$no_split,
             'help'       => \$help, );

if ( $help or (! $gene2cds) or (! $annot) ) { 
    die "Format: cds2gene_annot.pl\n",
        "        --gene2cds|-g  [gene-to-CDS table]\n",
        "        --annot|-a     [cds-annots in TSV]\n",
        "        --title|-t     [optional title of annot column (default, \"Annot\")]\n",
        "        --no_split|-n  [option to not split input annotations]\n",
        "        --help|-h      [print this message]\n",
        ;
}

if (! $title) {
    $title = 'Annot';
}

open my $GENE, '<', $gene2cds;
while (my $input = <$GENE>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $gene = $1;
        my $cds  = $2;
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else { 
        warn "From transcript-to-gene table $gene2cds, can't parse input: $input\n";
    }
}
close $GENE;

open my $ANNOT, '<', $annot;
while (my $input = <$ANNOT>) {
    chomp $input;
    if ( $input !~ /\A [#] /xms ) { 
        if ( $input =~ /\A (\S+) \t (.+) \z/xms ) { 
            my $cds         = $1;
            my $orig_annots = $2;

            if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
                die "Cannot map CDS $cds to a gene in: $input\n";
            }

            my $gene   = $data_ref->{'cds'}->{$cds}->{'gene'};

            $orig_annots =~ s/\A\s+//;
            $orig_annots =~ s/\s+\z//;

            if (! $no_split ) {
                my @annots = split /[,]|\t/, $orig_annots;

                foreach my $annot (@annots) {
                    $data_ref->{'gene'}->{$gene}->{'annots'}->{$annot} = 1;
                }
            }
            else {
                $data_ref->{'gene'}->{$gene}->{'annots'}->{$orig_annots} = 1;
            }
        }
        else { 
            warn "From annotation file $annot, can't parse input: $input\n";
        } 
    }
}
close $ANNOT;

my @genes = sort keys %{ $data_ref->{'gene'} };

my $header = "Gene\t$title\n";

foreach my $gene1 (@genes) {
    my $annot_text  = q{};
    my @annot_descs = ();
    if ( exists $data_ref->{'gene'}->{$gene1}->{'annots'} ) { 
        @annot_descs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'annots'} };
        $annot_text  = join '; ', @annot_descs;
        print $header if $header;
        $header = 0;
        print "$gene1\t$annot_text\n";
    }
}

