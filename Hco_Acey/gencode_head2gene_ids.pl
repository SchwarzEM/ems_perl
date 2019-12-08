#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tGene_name\tGene_type";

my $data_ref;

while (my $input = <>) {
    chomp $input;

    # Sample desired input lines; get this with "grep '>'".
    # >ENSMUST00000017290.10|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000005869.1|Brca1-201|Brca1|6572|protein_coding|
    # >ENSMUST00000156843.1|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000005871.1|Brca1-204|Brca1|2136|retained_intron|
    # >ENSMUST00000131460.1|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000005872.1|Brca1-202|Brca1|539|retained_intron|
    # >ENSMUST00000191198.1|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000119752.2|Brca1-207|Brca1|531|protein_coding|
    # >ENSMUST00000188168.1|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000119753.1|Brca1-205|Brca1|852|lncRNA|
    # >ENSMUST00000142086.2|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000005870.2|Brca1-203|Brca1|1917|nonsense_mediated_decay|
    # >ENSMUST00000190862.1|ENSMUSG00000017146.12|OTTMUSG00000002870.3|OTTMUST00000119754.2|Brca1-206|Brca1|355|protein_coding|

    if ( $input =~ /\A [>] 
                       (?:[^\s\|]+) \|  # transcript name
                       ([^\s\|]+)   \|  # Ensembl (Gencode) gene ID
                       (?:[^\s\|]+) \|  # OTTMUSG-type gene ID
                       (?:[^\s\|]+) \|  # OTTMUSG-type transcript ID 
                       (?:[^\s\|]+) \|  # human-readable transcript name
      	       	       ([^\s\|]+)   \| 	# Human-readable gene name
                       (?:\d+)      \|  # Mysterious integer -- transcript length in nt, perhaps?
                       ([^\s\|]+)   \|  # *Transcript* type (note that 2+ transcripts can merit 2+ types for one gene)
                    \z/xms ) {
        my $gene_id   = $1;
        my $gene_name = $2;
        my $tx_type   = $3;
        $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'}->{$gene_name} = 1;
        $data_ref->{'gene_id'}->{$gene_id}->{'tx_type'}->{$tx_type}     = 1;
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

my @gene_ids = sort keys %{ $data_ref->{'gene_id'} };

foreach my $gene_id (@gene_ids) {
    my @gene_names     = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'} };
    my $gene_name_list = join '|', @gene_names;

    my @tx_types     = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'tx_type'} };
    my $tx_type_list = join '|', @tx_types;

    # Print header once at top.
    print "$header\n" if $header;
    $header = q{};

    print "$gene_id\t$gene_name_list\t$tx_type_list\n";
}

