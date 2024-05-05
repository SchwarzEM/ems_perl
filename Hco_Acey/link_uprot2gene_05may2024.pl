#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gene2gb_data  = q{};
my $uprot2gb_data = q{};

my $data_ref;

my $header = "Gene\tGenBank_ID\tUniProt_IDs";

$gene2gb_data  = $ARGV[0] if $ARGV[0];
$uprot2gb_data = $ARGV[1] if $ARGV[1];

if ( (! $gene2gb_data ) and (! $uprot2gb_data ) ) {
    die "Format: link_uprot2gene_05may2024.pl",
        " [gene to GenBank ID table] [UniProt to GenBank ID table]",
        " > [gene ID / GenBank ID / UniProt IDs table]",
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

        # Mappings of gene names to GenBank IDs should be unique / 1-to-1:
        if ( exists $data_ref->{'gene'}->{$gene}->{'gb_id'} ) {
            die "Redundant mapping of gene $gene to GenBank IDs $data_ref->{'gene'}->{$gene}->{'gb_id'} and $gb_id\n";
        }
        else {
            $data_ref->{'gene'}->{$gene}->{'gb_id'} = $gb_id;
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

        # Note that a GenBank *gene* ID can map to two more more UniProt *protein* IDs.
        $data_ref->{'gb_id'}->{$gb_id}->{'uprot_id'}->{$uprot_id} = 1;
    }
    else {
        die "From UniProt ID to GenBank ID table $uprot2gb_data, cannot parse: $input\n";
    }
}
close $UPROT2GB_DATA;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    if (! exists $data_ref->{'gene'}->{$gene}->{'gb_id'} ) {
        die "Cannot map gene $gene to a GenBank ID\n";
    }
    my $gb_id = $data_ref->{'gene'}->{$gene}->{'gb_id'};

    # There exist a few instances of genes with GenBank IDs which, for whatever reason, were not given UniProt IDs;
    #    the entire mapping should not be killed because of this, but, such exceptions should be noted.
    if (! exists $data_ref->{'gb_id'}->{$gb_id}->{'uprot_id'} ) {
        warn "Cannot map gene $gene and GenBank ID $gb_id to UniProt IDs\n";
    }
    else {
        print "$header\n" if $header;
        $header = q{};

        my @uniprots = sort keys %{ $data_ref->{'gb_id'}->{$gb_id}->{'uprot_id'} };
        my $uprot_text = join '; ', @uniprots;
        print "$gene\t$gb_id\t$uprot_text\n";
    }
}

