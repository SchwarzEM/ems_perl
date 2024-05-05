#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gene2gb_data  = q{};
my $uprot2gb_data = q{};

my $data_ref;

my $header = "UniProt_ID\tGenBank_ID\tGene";

$gene2gb_data  = $ARGV[0] if $ARGV[0];
$uprot2gb_data = $ARGV[1] if $ARGV[1];

if ( (! $gene2gb_data ) and (! $uprot2gb_data ) ) {
    die "Format: link_uprot2gene_05may2024.pl",
        " [gene to GenBank ID table] [UniProt to GenBank ID table]",
        " > [ UniProt ID / GenBank ID / gene ID table]",
        "\n",
        ;
}

open my $GENE2GB_DATA, '<', $gene2gb_data;
while (my $input = <$GENE2GB_DATA>) {
    chomp $input;

    # Sample input:
    # Sherm_I.g354    QR680_007293
    # Sherm_I.g355    QR680_007294
    # Sherm_I.g356    QR680_007295

    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene  = $1;
        my $gb_id = $2;

        # Mappings of GenBank IDs to gene names should be unique / 1-to-1:
        if ( exists $data_ref->{'gb_id'}->{$gb_id}->{'gene'} ) {
            die "Redundant mapping of GenBank ID $gb_id to genes $data_ref->{'gb_id'}->{$gb_id}->{'gene'} and $gene\n";
        }
        else {
            $data_ref->{'gb_id'}->{$gb_id}->{'gene'} = $gene;
        }
    }
    else {
        die "From gene to GenBank ID table $gene2gb_data, cannot parse: $input\n";
    }
}
close $GENE2GB_DATA;

open my $UPROT2GB_DATA, '<', $uprot2gb_data;
while (my $input = <$UPROT2GB_DATA>) {
    chomp $input;

    # Sample input:
    # A0AA39H2E1      QR680_002165
    # A0AA39H3B4      QR680_002165
    # A0AA39HQU0      QR680_004766

    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $uprot_id = $1;
        my $gb_id    = $2;

        # Note that UniProt protein IDs *are* allowed to map to 2+ GenBank IDs, because a 2+ genes can encode identical peptides!
        $data_ref->{'uprot_id'}->{$uprot_id}->{'gb_id'}->{$gb_id} = 1;
    }
    else {
        die "From UniProt ID to GenBank ID table $uprot2gb_data, cannot parse: $input\n";
    }
}
close $UPROT2GB_DATA;

my @uniprots = sort keys %{ $data_ref->{'uprot_id'} };
foreach my $uprot_id (@uniprots) {
    if (! exists  $data_ref->{'uprot_id'}->{$uprot_id}->{'gb_id'} ) {
        die "Cannot map UniProt ID $uprot_id to a GenBank ID\n";
    }
    my @gb_ids = sort keys %{ $data_ref->{'uprot_id'}->{$uprot_id}->{'gb_id'} };
    foreach my $gb_id (@gb_ids) {
        if (! exists $data_ref->{'gb_id'}->{$gb_id}->{'gene'} ) {
            die "Cannot map UniProt ID $uprot_id and GenBank ID $gb_id to original gene name\n";
        }
        else {
            print "$header\n" if $header;
            $header = q{};

            my $gene = $data_ref->{'gb_id'}->{$gb_id}->{'gene'};
            print "$uprot_id\t$gb_id\t$gene\n";
        }
    }
}

