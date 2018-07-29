#!/usr/bin/env perl

use strict;
use warnings;
use autodie ;

my $infile = q{};
my $i      = 1;

$infile = $ARGV[0] if $ARGV[0];
$i      = $ARGV[1] if $ARGV[1];

my $data_ref;

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t \S+ \t \S+ \t \S+ \t ([^\t]+) \t /xms ) { 
        my $gene = $1;
        my $phen = $2;
        if ( $phen =~ /\A \S+ \s+ lethal \z/xms ) {
            $data_ref->{'let_gene'}->{$gene}->{'let_phen'}->{$phen} = 1;
        }
        $data_ref->{'gene'}->{$gene}->{'annot'}->{$input} = 1;
    }
}
close $INFILE;

my @genes = sort %{ $data_ref->{'let_gene'} };
foreach my $gene (@genes) {
     my $let_count = keys %{ $data_ref->{'let_gene'}->{$gene}->{'let_phen'} };
     if ( $let_count >= $i ) {
         my @annots = sort keys %{ $data_ref->{'gene'}->{$gene}->{'annot'} };
         foreach my $annot (@annots) {
             print "$annot\n";
         }
    }
}

