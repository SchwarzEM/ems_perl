#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use JSON::PP;

my @inputs = @ARGV;

my $data_ref;

if (! @inputs ) {
    die "Format: afserv_json_22nov2025.pl [job_name] [1+ FASTAs] > job_name.json ;\n"; 
}

# Use this array to keep the order of input sequences unchanged in the JSON.
my @input_seqs = ();

my $job = shift @inputs;

foreach my $fasta (@inputs) {
    my $seq = q{};
    open my $FASTA, '<', $fasta;
    while ( my $input = <$FASTA> ) {
        chomp $input;
        if ( $input =~ /\A [>] (\S+) /xms ) {
            $seq = $1;
            push @input_seqs, $seq;
            $data_ref->{'seq'}->{$seq}->{'res'} = q{};
        }
        elsif ( $input =~ /\S/xms ) {
            $input =~ s/\s//g;
            $data_ref->{'seq'}->{$seq}->{'res'} .= $input;
        }
    }
    close $FASTA;
}

my @sequence_entries = map {
    {
        proteinChain => {
            sequence => $data_ref->{'seq'}->{$_}->{'res'},
            count    => 1,
        }
    }
} @input_seqs;

my $json_data = [
    {
        name        => $job,
        modelSeeds  => [],
        sequences   => \@sequence_entries,
        dialect     => "alphafoldserver",
        version     => 1,
    }
];

# Rename this for file naming rather than internal job naming.
$job    = "$job.json";
$job    = safename($job);

open my $JSON_FILE, '>', $job;
my $json = JSON::PP->new->utf8->canonical->pretty->encode($json_data);
print $JSON_FILE $json;
close $JSON_FILE;

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

