#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tGene_Ontology";

my $go_terms    = q{};
my $ens2mgi_ids = q{};
my $mgi2go_ids  = q{};

my $id            = q{};
my $desc          = q{};

my %go_id2go_desc = ();

my %mgi2ens       = ();

$go_terms    = $ARGV[0] if $ARGV[0];
$ens2mgi_ids = $ARGV[1] if $ARGV[1];
$mgi2go_ids  = $ARGV[2] if $ARGV[2];

open my $GO_TERMS, '<', $go_terms;
while ( my $input = <$GO_TERMS> ) {
    chomp $input;
    if ( $input =~ /\A id [:] [ ] (GO[:]\d+) /xms ) {
        $id = $1;
    }
    elsif ( ( $input =~ /\A name[:] [ ] (\S.+\S) \s* \z/xms ) and ($id) ) {
        $desc = $1;

        $go_id2go_desc{$id} = $desc;

        $id   = q{};
        $desc = q{};
    }
}
close $GO_TERMS;

open my $ENS2MGI, '<', $ens2mgi_ids;
while ( my $input = <$ENS2MGI> ) {
    chomp $input;
    if ( $input =~ /\A (ENSMUSG\d+) \t (MGI[:]\d+) \z/xms ) {
        my $ens_id = $1;
        my $mgi_id = $2;
        $mgi2ens{$mgi_id} = $ens_id;
    }
    else {
       die "From $ens2mgi_ids, cannot parse: $input\n";
    }
}

close $ENS2MGI;

open my $MGI2GO, '<', $mgi2go_ids;
while ( my $input = <$MGI2GO> ) {
    chomp $input;
    if ( $input =~ /\A (MGI[:]\d+) \t (GO[:]\d+) \z/xms ) {
        my $mgi_id = $1;
        my $go_id  = $2;
        if ( ( exists $mgi2ens{$mgi_id} ) and ( exists $go_id2go_desc{$go_id} ) ) {
            my $gene    = $mgi2ens{$mgi_id};
            my $go_desc = $go_id2go_desc{$go_id};
            $go_desc    = "$go_desc [$go_id]";
            $data_ref->{'gene'}->{$gene}->{'go_desc'}->{$go_desc} = 1;
        }
    }
    else {
       die "From $mgi2go_ids, cannot parse: $input\n";
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @go_descs = sort keys %{ $data_ref->{'gene'}->{$gene}->{'go_desc'} };
    my $go_desc_text = join '; ', @go_descs;

    print "$header\n" if $header;
    $header = q{};

    print "$gene\t$go_desc_text\n";
}

