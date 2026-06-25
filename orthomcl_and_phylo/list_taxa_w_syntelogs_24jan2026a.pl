#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $infile    = q{};
my @pangenes  = ();

my $data_ref;

$infile    = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: list_taxa_w_syntelogs_24jan2026a.pl [infile] > [taxon set; list of pangenes with that taxon set]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;

    # Sample input:
    # ofID    pgID    interpChr       interpOrd       pgRepID genome  og      flag    id      chr     start   end     ord
    # 0_368   1       Necator_chrI    1       0_368   Aroian  35669   PASS    Necator_chrI.g137       Necator_chrI    13892   14485   1

    if ( ( $input !~ /\A ofID \t/xms ) and 
         ( $input =~ /\A \S+ \t (\S+) \t (?: \S+ \t){3} (\S+) \t (?: \S+ \t) (\S+) \t \S+ \t .* \z/xms ) ) {
        my $pangene = $1;
        my $taxon   = $2;
        my $flag    = $3;
        $data_ref->{'pangene'}->{$pangene}->{'flag'}->{$flag}->{'taxon'}->{$taxon} = 1;
        push @pangenes, $pangene;
    }
}
close $INFILE;

@pangenes = uniq(@pangenes);

foreach my $pangene (@pangenes) {
    if ( exists $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'} ) {
        my @taxa = sort keys %{ $data_ref->{'pangene'}->{$pangene}->{'flag'}->{'PASS'}->{'taxon'} };
        my $taxa_text = join '; ', @taxa;
        $data_ref->{'taxa_text'}->{$taxa_text}->{'pangene'}->{$pangene} = 1;
    }
}

my @taxa_texts = sort keys %{ $data_ref->{'taxa_text'} };

foreach my $taxa_text (@taxa_texts) {
    my @pangenes = sort keys %{ $data_ref->{'taxa_text'}->{$taxa_text}->{'pangene'} };
    my $pangene_text =  join '; ', @pangenes;
    my $pangene_count = @pangenes;
    print "$taxa_text\t$pangene_count\t$pangene_text\n";
}

