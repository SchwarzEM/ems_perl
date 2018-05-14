#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile = $ARGV[0];
my $count  = $ARGV[1];

my $data_ref;

my $header = q{};

open my $INFILE, '<', $infile;
while (my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A Gene /xms ) { 
        $header = $input;
    }
    elsif ( $input =~ /\A (\S+) \t /xms ) { 
        my $gene = $1;
        $data_ref->{'gene'}->{$gene}->{'count'}++;
        $data_ref->{'gene'}->{$gene}->{'annot'}->{$input} = 1;
    }
    else { 
        die "Cannot parse input: $input\n";
    }
}
close $INFILE;

my @genes = sort keys %{ $data_ref->{'gene'} };
@genes    = sort { $data_ref->{'gene'}->{$b}->{'count'} <=> $data_ref->{'gene'}->{$a}->{'count'} } @genes;

$count--;

my @select_genes = @genes[0..$count];
@select_genes    = sort @select_genes;

foreach my $gene (@select_genes) {
    print $header if $header;
    $header = q{};
    my @annots = sort keys %{  $data_ref->{'gene'}->{$gene}->{'annot'} };
    @annots = sort @annots;
    @annots = uniq(@annots);
    foreach my $annot (@annots) {
        print "$annot\n";
    }
}

