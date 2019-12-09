#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Gene\tMutant_phenotypes";

my $pheno_terms  = q{};
my $pheno_annots = q{};

my $id      = q{};
my $desc    = q{};
my %id2desc = ();

$pheno_terms  = $ARGV[0] if $ARGV[0];
$pheno_annots = $ARGV[1] if $ARGV[1];

open my $PHENO_TERMS, '<', $pheno_terms;
while ( my $input = <$PHENO_TERMS> ) {
    chomp $input;
    if ( $input =~ /\A id [:] [ ] (MP[:]\d+) /xms ) {
        $id = $1;
    }
    elsif ( ( $input =~ /\A name[:] [ ] (\S.+\S) \s* \z/xms ) and ($id) ) {
        $desc = $1;

        $id2desc{$id} = $desc;

        $id   = q{};
        $desc = q{};
    }
}
close $PHENO_TERMS;

open my $PHENO_ANNOTS, '<', $pheno_annots;

while ( my $input = <$PHENO_ANNOTS> ) {
    chomp $input;
    if ( $input	=~ /\A (ENSMUSG\d+) \t (\S+) \z/xms ) {
       	my $id            = $1;
        my $pheno_acc_txt = $2;

        my @pheno_accs      = split /,/, $pheno_acc_txt;
        my @pheno_annots    = ();
        my $pheno_annot_txt = q{};

        foreach my $pheno_acc (@pheno_accs) {
            my $pheno_desc = $id2desc{$pheno_acc};
            $pheno_desc = "$pheno_desc [$pheno_acc]";
            push @pheno_annots, $pheno_desc;
        }
        $pheno_annot_txt = join '; ', @pheno_annots;

        print "$header\n" if $header;
        $header = q{};

        print "$id\t$pheno_annot_txt\n";
    }
    else {
        die "From $pheno_annots, cannot parse: $input\n";
    }
}
close $PHENO_ANNOTS;




my @ens_ids = sort keys %{ $data_ref->{'ens_id'} };
foreach my $ens_id (@ens_ids) {
    my $symbol = $data_ref->{'ens_id'}->{$ens_id}->{'symbol'};
    my $desc   = $data_ref->{'ens_id'}->{$ens_id}->{'desc'};

    print "$header\n" if $header;
    $header = q{};
    print "$ens_id\t$symbol\t$desc\n";
}


