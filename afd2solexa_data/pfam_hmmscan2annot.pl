#!/usr/bin/env perl

# pfam_hmmscan2annot.pl -- Erich Schwarz <ems394@cornell.edu>, orig. 9/2/2016; sig. revised 2/18/2026.
# Purpose: given either a gene2cds or cds2gene* table, and a PFAM/hmmscan --tblout or --domtblout output, make a 2-column table of genes and PFAM domains.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $cds2gene  = q{};
my $gene2cds  = q{};
my $pfam_gene = q{};
my $pfam_dom  = q{};
my $help;

GetOptions ( 'cds2gene=s'  => \$cds2gene,
             'gene2cds=s'  => \$gene2cds,
             'pfam_gene=s' => \$pfam_gene,
             'pfam_dom=s'  => \$pfam_dom,
             'help'        => \$help, );

if (    $help 
     or ( (! $cds2gene ) and (! $gene2cds ) ) 
     or ( $cds2gene and $gene2cds ) 
     or ( (! $pfam_gene ) and (! $pfam_dom ) ) 
     or ( $pfam_gene and $pfam_dom )
   ) { 
    die "Format: pfam_hmmscan2annot.pl\n",
        "            --cds2gene|-c [CDS-to-gene table] or --gene2cds|-g [gene-to-CDS table]\n",
        "            --pfam_gene [PFAM/hmmscan --tblout tabular output] or --pfam_dom [PFAM/hmmscan --domtblout tabular output]\n",
        "            --help|-h [print this message]\n",
        ;
}

my $GENE_CDS;
if ($cds2gene) {
    open $GENE_CDS, '<', $cds2gene;
}
elsif ($gene2cds) {
    open $GENE_CDS, '<', $gene2cds;
}
while ( my $input = <$GENE_CDS> ) { 
    chomp $input;
    my $cds  = q{};
    my $gene = q{};
    if ( $cds2gene and ( $input =~ /\A (\S+) \t (\S+) \z/xms ) ) { 
        $cds  = $1;
        $gene = $2;
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    elsif ( $gene2cds and ( $input =~ /\A (\S+) \t (\S+) \z/xms ) ) {
        $gene = $1;
        $cds  = $2;
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else { 
        die "From transcript-to-gene table \"$cds2gene\" or gene-to-transcript table \"$gene2cds\", can't parse input: $input\n";
    }
}
close $GENE_CDS;

my $pfam_input = q{};

if ( $pfam_gene ) {
    $pfam_input = $pfam_gene;
}
elsif ( $pfam_dom ) {
    $pfam_input = $pfam_dom;
}
else {
    die "Cannot identify either Pfam/gene ($pfam_gene) or Pfam/domain ($pfam_dom) table as input.\n";
}

open my $PFAM, '<', $pfam_input;
while (my $input = <$PFAM>) {
    chomp $input;
    if ( $input !~ /\A [#] /xms ) { 
        my $pfam_name = q{};
        my $pfam_acc  = q{};
        my $cds       = q{};

        # Sample input line from pfam_gene:
        # HlyIII    PF03006.27    Necator_2022.05.29.01.07.g2.t3    -  [...]
        if ( ( $pfam_gene ) and ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ /xms ) ) { 
            $pfam_name = $1;
            $pfam_acc  = $2;
            $cds       = $3;
            if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
                die "Cannot map CDS $cds to a gene\n";
            }
            my $gene      = $data_ref->{'cds'}->{$cds}->{'gene'};
            my $pfam_desc = "$pfam_name [$pfam_acc]";
            $data_ref->{'gene'}->{$gene}->{'pfam_desc'}->{$pfam_desc} = 1;
        }
        # Sample input line from pfam_dom: 
        # HlyIII    PF03006.27    224    Necator_2022.05.29.01.07.g2.t3    -  [...]
        elsif ( ( $pfam_dom ) and ( $input =~ /\A (\S+) \s+ (\S+) \s+ \d+ \s+ (\S+) \s+ /xms ) ) {
            $pfam_name = $1;
            $pfam_acc  = $2;
            $cds       = $3;
            if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
                die "Cannot map CDS $cds to a gene\n";
            }
            my $gene      = $data_ref->{'cds'}->{$cds}->{'gene'};
            my $pfam_desc = "$pfam_name [$pfam_acc]";
            $data_ref->{'gene'}->{$gene}->{'pfam_desc'}->{$pfam_desc} = 1;
        }
        else { 
            die "From PFAM table $pfam_input, can't parse input: $input\n";
        } 
    }
}
close $PFAM;

my @genes = sort keys %{ $data_ref->{'gene'} };

my $header = "Gene\tPFAM\n";

foreach my $gene1 (@genes) {
    my $pfam_text  = q{};
    my @pfam_descs = ();
    if ( exists $data_ref->{'gene'}->{$gene1}->{'pfam_desc'} ) { 
        @pfam_descs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'pfam_desc'} };
        $pfam_text  = join '; ', @pfam_descs;
    }
    print $header if $header;
    $header = 0;
    print "$gene1\t$pfam_text\n";
}

