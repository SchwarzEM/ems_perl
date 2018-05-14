#!/usr/bin/env perl

my $gene_index  = $ARGV[0];
my $olist_index = $ARGV[1];
my $omim_index  = $ARGV[2];

my %ensg2hgnc = ();
my %ortholist = ();
my %omimlist  = ();

open my $GENES, '<', $gene_index or die "Can't open gene index $gene_index: $!";
while (my $input = <$GENES>) {
    chomp $input;
    # Sample input:
    # human|ENSG00000198712|MT-CO2
    if ( $input =~ /\A human \| (ENSG\d+) \| (\S+) \z/xms ) { 
        my $ensg          = $1;
        my $hgnc          = $2;
        $ensg2hgnc{$ensg} = $hgnc;
    }
    else { 
        warn "From gene index $gene_index, can't parse input line: $input\n";
    }
}
close $GENES or die "Can't close filehandle to gene index $gene_index: $!";

open my $OLIST, '<', $olist_index or die "Can't open OrthoList human gene list $olist_index: $!";
while (my $input = <$OLIST>) {
    chomp $input;
    # Sample input:
    # ENSG00000000003
    if ( $input =~ /\A (ENSG\d{11}) \z/xms ) {
        my $ensg          = $1;
        $ortholist{$ensg} = 1;
    }
    else {
        die "From OrthoList human gene list $olist_index, can't parse input line: $input\n";
    }
}
close $OLIST or die "Can't close filehandle to OrthoList human gene list $olist_index: $!";

open my $OMIM, '<', $omim_index or die "Can't open OMIM morbid cloned human disease gene list $omim_index: $!";
while (my $input = <$OMIM>) {
    chomp $input;
    # Sample input:
    # 1C7
    # 3H11AG
    # 3M2
    if ( $input =~ /\A (\S+) \s* \z/xms ) {
        my $hgnc = $1;
        $omimlist{$hgnc} = 1;
    }
    else {
        die "From OMIM morbid cloned human disease gene list $omim_index, can't parse input line: $input\n";
    }
}
close $OMIM or die "Can't close filehandle to OMIM morbid cloned human disease gene list $omim_index: $!";

my @olist_hgncs = map { $ensg2hgnc{$_} } sort keys %ortholist;
my @omim_hgncs  = sort keys %omimlist;
my @overlaps    = ();

foreach my $olist_hgnc (@olist_hgncs) { 
    if ( exists $omimlist{$olist_hgnc} ) { 
        push @overlaps, $olist_hgnc;
    }
}

my $olist_hgnc_count = @olist_hgncs;
my $omim_hgnc_count  = @omim_hgncs;
my $overlap_count    = @overlaps;

print "OrthoList HGNCs:    $olist_hgnc_count\n";
print "OMIM HGNCs:         $omim_hgnc_count\n";
print "Overlapping HGNCs:  $overlap_count\n";

