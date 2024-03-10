#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $proteome    = q{};
$proteome       = $ARGV[0] if $ARGV[0];

my $prot_id     = q{};
my $prot_seq    = q{};

if (! $proteome ) {
    die "Format: prots_to_nonred_seqs_11jul2023.pl [proteome] => [seq name (tab) aa seq in one text line]\n";
}

open my $PROTEOME, '<', $proteome;
while (my $input = <$PROTEOME> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) {
        my $new_prot_id = $1;
        if ( $prot_id and $prot_seq ) {
            print "$prot_id\t$prot_seq\n";
        }
        $prot_id  = $new_prot_id;
        $prot_seq = q{};
    }
    elsif ($prot_id) {
        $input =~ s/\s//;
        $prot_seq = $prot_seq . $input;
    }
}
close $PROTEOME;

# Finish recording any remaining sequence data:
if ( $prot_id and $prot_seq ) {
    print "$prot_id\t$prot_seq\n";
}

