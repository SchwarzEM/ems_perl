#!/usr/bin/env perl

# censor_gff_nofastas.pl -- Erich Schwarz <ems@emstech.org>, 2/9/2014.
# Purpose: given a FASTA sequence file and a GFF, censor any non-comment line that invokes a sequence name not actually present in the FASTA.

use strict;
use warnings;

my $fasta = $ARGV[0];
my $gff   = $ARGV[1];

my %obs = ();

if ( (! $fasta) or (! $gff) ) {
    die "Format: censor_gff_nofastas.pl [FASTA] [GFF] > [GFF censored of non-comment lines that invoke sequence names not actually present in the FASTA] ;\n";
}

open my $FASTA, '<', $fasta or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $seq = $1;
        $obs{$seq} = 1;
    }
}
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

open my $GFF, '<', $gff or die "Can't open GFF file $gff: $!";
while (my $input = <$GFF>) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A (\S+) /xms ) ) {
        my $seq = $1;
        print "$input\n" if (exists $obs{$seq});
    }
    else { 
        print "$input\n";
    }
}
close $GFF or die "Can't close filehandle to GFF file $gff: $!";

