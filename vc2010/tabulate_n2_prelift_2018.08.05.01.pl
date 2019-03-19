#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Gene\tPre_liftover\n";

my $well_trans  = q{};
my $badly_trans = q{};

my $data_ref;

$well_trans  = $ARGV[0] if $ARGV[0];
$badly_trans = $ARGV[1]	if $ARGV[1];

open my $WELL, '<', $well_trans;
while (my $input = <$WELL>) {
    chomp $input;
    $data_ref->{'gene'}->{$input} = 'Translated';
}
close $WELL;

open my $BADLY, '<', $badly_trans;

while (my $input = <$BADLY>) {
    chomp $input;
    $data_ref->{'gene'}->{$input} = 'Mistranslated';
}
close $BADLY;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    print $header if $header;
    $header = q{};
    my $annot = $data_ref->{'gene'}->{$gene};
    print "$gene\t$annot\n";
}


