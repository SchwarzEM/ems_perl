#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $gene2cds = q{};
my $genelist = q{};

my $data_ref;

my @proteins = ();

$gene2cds = $ARGV[0] if $ARGV[0];
$genelist = $ARGV[1] if $ARGV[1];

if ( (! $gene2cds ) or (! $genelist ) ) {
    die "Format: flygene2prot_21aug2024.pl [gene2cds] [gene list] > [protein list]\n";
}

open my $GENE2CDS, '<', $gene2cds;
while ( my $input = <$GENE2CDS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene    = $1;
        my $protein = $2;
        $data_ref->{gene}->{$gene}->{protein}->{$protein} = 1;
    }
    else {
        die "From gene2cds file $gene2cds, cannot parse: $input\n";
    }
}
close $GENE2CDS;

open my $GENELIST, '<', $genelist;
while ( my $input = <$GENELIST> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $gene = $1;
        if ( exists $data_ref->{gene}->{$gene}->{protein} ) { 
            my @proteins1 = sort keys %{ $data_ref->{gene}->{$gene}->{protein} };
            push @proteins, @proteins1;
        }
    }
    else {
        die "From genelist file $genelist, cannot parse: $input\n";
    }
}
close $GENELIST;

@proteins = sort(@proteins);
@proteins = uniq(@proteins);

foreach my $protein (@proteins) {
    print "$protein\n";
}
