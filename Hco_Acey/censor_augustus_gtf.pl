#!/usr/bin/env perl

use strict;
use warnings;

my $bad_genes = $ARGV[0];
my $bad_txs   = $ARGV[1];
my $orig_aug  = $ARGV[2];

my $data_ref;

if ( (! $bad_genes) or (! $bad_txs) or (! $orig_aug) ) {
    die "Format: censor_augustus_gtf.pl [bad gene list] [bad transcript list] [original AUGUSTUS GTF2] > [new AUGUSTUS GTF2]\n";
}

open my $BAD_GENES, '<', $bad_genes or die "Can't open bad gene list $bad_genes: $!";
while (my $input = <$BAD_GENES>) { 
    chomp $input;
    $data_ref->{'bad_gene'}->{$input} = 1;
}
close $BAD_GENES or die "Can't close filehandle to bad gene list $bad_genes: $!";

open my $BAD_TXS, '<', $bad_txs or die "Can't open bad transcript list $bad_txs: $!";
while (my $input = <$BAD_TXS>) {
    chomp $input;
    $data_ref->{'bad_tx'}->{$input} = 1;
}
close $BAD_TXS or die "Can't close filehandle to bad transcript list $bad_txs: $!";

open my $ORIG_AUG, '<', $orig_aug or die "Can't open original AUGUSTUS GTF2 file $orig_aug: $!";

while (my $input = <$ORIG_AUG>) { 
    chomp $input;
    my $print = 1;

    # Sample input lines:
    # Acey_s0001_scaf AUGUSTUS        gene    1       7875    0.3     -       .       Acey_s0001.g1
    # Acey_s0001_scaf AUGUSTUS        transcript      1       7875    0.3     -       .       Acey_s0001.g1.t1
    # Acey_s0001_scaf AUGUSTUS        intron  1       1284    0.3     -       .       transcript_id "Acey_s0001.g1.t1"; gene_id "Acey_s0001.g1";

    # Sieve out lines that are disqualified in any of three ways.
    if ( ( $input =~ /\A (?: [^\t]* \t){2} gene \t (?: [^\t]* \t){5} (\S+) /xms ) and ( exists $data_ref->{'bad_gene'}->{$1} ) ) { 
        $print = 0;
    }
    if ( ( $input =~ /\A (?: [^\t]* \t){2} transcript \t (?: [^\t]* \t){5} (\S+) /xms ) and ( exists $data_ref->{'bad_tx'}->{$1} ) ) {
        $print = 0;
    }  
    if ( ( $input =~ / transcript_id [ ] \" ( [^\s\"]+ ) \" ; /xms ) and ( exists $data_ref->{'bad_tx'}->{$1} ) ) {
        $print = 0;
    }

    # Print any lines that survive the sieving.
    if ($print) { 
        print "$input\n";
    }
}
close $ORIG_AUG or die "Can't close filehandle to original AUGUSTUS GTF2 file $orig_aug: $!";

