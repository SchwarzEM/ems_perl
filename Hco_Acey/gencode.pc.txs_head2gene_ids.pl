#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tGene_name";

my $data_ref;

while (my $input = <>) {
    chomp $input;

    # Sample desired input lines; get this with "grep '>'".
    # >ENSMUST00000070533.5|ENSMUSG00000051951.6|OTTMUSG00000026353.2|OTTMUST00000065166.1|Xkr4-201|Xkr4|3634|UTR5:1-150|CDS:151-2094|UTR3:2095-3634|
    # >ENSMUST00000208660.2|ENSMUSG00000025900.14|OTTMUSG00000049985.3|OTTMUST00000145515.1|Rp1-202|Rp1|4170|UTR5:1-54|CDS:55-4170|
    # >ENSMUST00000027032.6|ENSMUSG00000025900.14|OTTMUSG00000049985.3|OTTMUST00000127195.2|Rp1-201|Rp1|6869|UTR5:1-127|CDS:128-6415|UTR3:6416-6869|
    # >ENSMUST00000027035.10|ENSMUSG00000025902.14|OTTMUSG00000050014.7|OTTMUST00000127245.2|Sox17-201|Sox17|3127|UTR5:1-1082|CDS:1083-2342|UTR3:2343-3127|
    # >ENSMUST00000195555.2|ENSMUSG00000025902.14|OTTMUSG00000050014.7|OTTMUST00000127249.1|Sox17-208|Sox17|1977|UTR5:1-635|CDS:636-1511|UTR3:1512-1977|
    # >ENSMUST00000192650.6|ENSMUSG00000025902.14|OTTMUSG00000050014.7|OTTMUST00000127247.2|Sox17-206|Sox17|3242|UTR5:1-1851|CDS:1852-2916|UTR3:2917-3242|

    if ( $input =~ /\A [>] 
                       (?:[^\s\|]+) \|  # transcript name
                       ([^\s\|]+)   \|  # Ensembl (Gencode) gene ID
                       (?:[^\s\|]+) \|  # OTTMUSG-type gene ID
                       (?:[^\s\|]+) \|  # OTTMUSG-type transcript ID 
                       (?:[^\s\|]+) \|  # human-readable transcript name
      	       	       ([^\s\|]+)   \| 	# Human-readable gene name
                       (?:\d+)      \|  # Mysterious integer -- transcript length in nt, perhaps?
                    .* \z/xms ) {
        my $gene_id   = $1;
        my $gene_name = $2;
        $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'}->{$gene_name} = 1;
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

my @gene_ids = sort keys %{ $data_ref->{'gene_id'} };

foreach my $gene_id (@gene_ids) {
    my @gene_names     = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'gene_name'} };
    my $gene_name_list = join '|', @gene_names;

    # Print header once at top.
    print "$header\n" if $header;
    $header = q{};

    print "$gene_id\t$gene_name_list\n";
}

