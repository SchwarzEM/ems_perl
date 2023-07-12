#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $proteome    = q{};
$proteome       = $ARGV[0] if $ARGV[0];

my $prot_table  = q{};

my $prot_id     = q{};
my $prot_seq    = q{};

my $data_ref;

if (! $proteome ) {
    die "Format: prots_to_nonred_seqs_11jul2023.pl [proteome] => [lines of identical proteins]\n";
}

$prot_table = "$proteome.prots.txt";
$prot_table = safename($prot_table);

open my $PROTEOME, '<', $proteome;
while (my $input = <$PROTEOME> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) {
        my $new_prot_id = $1;
        if ( $prot_id and $prot_seq ) {
            $data_ref->{'prot_seq'}->{$prot_seq}->{'prot_id'}->{$prot_id} = 1;
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
    $data_ref->{'prot_seq'}->{$prot_seq}->{'prot_id'}->{$prot_id} = 1;
}

my @prot_seqs = keys %{ $data_ref->{'prot_seq'} };
foreach my $prot_seq1 (@prot_seqs) {
    my @prot_ids = sort keys %{ $data_ref->{'prot_seq'}->{$prot_seq1}->{'prot_id'} };
    my $prot_id_text = join '; ', @prot_ids;
    print "$prot_id_text\n";
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

