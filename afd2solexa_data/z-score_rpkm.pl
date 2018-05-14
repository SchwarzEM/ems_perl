#!/usr/bin/env perl

# z-score_rpkm.pl -- Erich Schwarz <emsch@caltech.edu>, 3/8/2012.
# Purpose: given a table with mean and sd, and an rpkm file, try to get z-scores of each gene in the rpkm file.

use strict;
use warnings;

my $statfile = $ARGV[0];
my $rpkmfile = $ARGV[1];
my $label    = $ARGV[2];

if (! ( $statfile and $rpkmfile and $label ) ) { 
    die "Usage: z-score_rpkm.pl [STATFILE] [RPKMFILE] [LABEL] > stdout ;\n";
}

my $wbgene  = q{};
my $mean    = q{};
my $sd      = q{};
my $rpkmval = q{};

my $data_ref;

open my $STAT, '<', $statfile or die "Can't open stat file $statfile: $!";
while (my $input = <$STAT>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S*) \s (\S+) \s (\S+) /xms ) { 
        $wbgene = $1;
        $mean   = $2;
        $sd     = $3;
        $data_ref->{'gene'}->{$wbgene}->{'mean'} = $mean;
        $data_ref->{'gene'}->{$wbgene}->{'sd'}   = $sd;
    }
}
close $STAT or die "Can't close filehandle to stat file $statfile: $!";

open my $RPKM, '<', $rpkmfile or die "Can't open rpkm file $rpkmfile: $!";
while (my $input = <$RPKM>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S*) \s \S+ \s (\S+) /xms ) { 
        $wbgene  = $1;
        $rpkmval = $2;
        if ( exists $data_ref->{'gene'}->{$wbgene} ) {
            $data_ref->{'gene'}->{$wbgene}->{'rpkmval'} = $rpkmval;
        }
    }
}
close $RPKM or die "Can't close filehandle to rpkm file $rpkmfile: $!";

print "Gene\tz-score_$label\n";

foreach my $stat_gene ( sort keys %{ $data_ref->{'gene'} } ) { 
    print $stat_gene;
    print "\t";
    my $z_score = 'n/a';
    if ( $data_ref->{'gene'}->{$stat_gene}->{'sd'} > 0 ) { 
        my $r_val = 0;
        if ( exists $data_ref->{'gene'}->{$stat_gene}->{'rpkmval'} ) { 
            $r_val = $data_ref->{'gene'}->{$stat_gene}->{'rpkmval'};
        }
        $z_score = ( ( $r_val - $data_ref->{'gene'}->{$stat_gene}->{'mean'} ) / $data_ref->{'gene'}->{$stat_gene}->{'sd'} );
        $z_score = sprintf("%.2f", $z_score);
    }
    print $z_score;
    print "\n";
}

