#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $genelist = q{};
my $proteome = q{};

$proteome = $ARGV[0] if $ARGV[0];
$genelist = $ARGV[1] if $ARGV[1];

if ( (! $proteome ) or (! $genelist ) ) {
    die "Format: xenprot_subset_21aug2024.pl [Xenopus proteome] [Xenopus gene list] > [Xenopus protein list]\n";
}

my $data_ref;

my @proteins = ();

open my $PROTEOME, '<', $proteome;
while ( my $input = <$PROTEOME> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) \s+ (\S+) \s* \z/xms ) {
        my $protein = $1;
        my $gene    = $2;
        $data_ref->{'gene'}->{$gene}->{'protein'}->{$protein} = 1;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "In proteome $proteome, cannot parse FASTA header line: $input\n";
    }
}
close $PROTEOME;

open my $GENELIST, '<', $genelist;
while ( my $input = <$GENELIST> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        my $gene    = $1;
        if ( exists $data_ref->{'gene'}->{$gene}->{'protein'} ) {
            my @proteins1 = sort keys %{ $data_ref->{'gene'}->{$gene}->{'protein'} };
            push @proteins, @proteins1;
        }
        else {
            warn "Cannot map gene $gene to a protein\n";
        }
    }
    else {
        die "In gene list $genelist, cannot parse: $input\n";
    }
}
close $GENELIST;

@proteins = sort(@proteins);
@proteins = uniq(@proteins);

foreach my $protein (@proteins) {
    print "$protein\n";
}
