#!/usr/bin/env perl

# tabulate_cegma_03nov2013.pl -- Erich Schwarz <ems394@cornell.edu>, 11/3/2013.
# Purpose: summarize several CEGMA reports in a row.  Requires slight hand-editing of the CEGMA completeness report texts.

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $data_ref;

my @stat_types = ( '#Prots', '%Completeness', '#Total',  'Average',  '%Ortho', );

my $dna  = q{};
my $type = q{};

my @dnas  = ();
my @types = qw(Complete Partial);

while (my $input = <>) {
    chomp $input;
    # Have to hand-edit lines like this in:
    # [Example:] DNA: C. elegans WS240 genome
    if ( $input =~ /\A DNA: \s+ (\S .+ \S) \s* \z/xms ) { 
        $dna = $1;
        push @dnas, $dna;
    }
    elsif ( $input =~ /\A \s+ (Complete|Partial) \s+ (\S .+ \S) \s* \z/xms ) {
        my $type       = $1;
        my $value_text = $2;
        my @values     = grep { $_ ne '-' } split /\s+/, $value_text;
        foreach my $i (0..4) {
            my $s_type = $stat_types[$i];
            $data_ref->{'dna'}->{$dna}->{$type}->{$s_type} = $values[$i];
        }
    }
}

@dnas = uniq @dnas;

print "Match\tValue";
foreach my $dna1 (@dnas) { 
    print "\t$dna1";
}
print "\n";

foreach my $type1 (@types) { 
    foreach my $s_type1 (@stat_types) { 
        print "$type1\t$s_type1";
        foreach my $dna2 (@dnas) {
            my $dna2_value = $data_ref->{'dna'}->{$dna2}->{$type1}->{$s_type1};
            print "\t$dna2_value";
        }
        print "\n";
    }
}

