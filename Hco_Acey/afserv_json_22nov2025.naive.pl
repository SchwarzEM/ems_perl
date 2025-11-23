#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @inputs = @ARGV;

my $data_ref;

if (! @inputs ) {
    die "Format: afserv_json_22nov2025.pl [job_name] [1+ FASTAs] > job_name.json ;\n"; 
}

my $job = shift @inputs;

my $header = q{};
$header    = q([) . "\n"
           . q(  {) . "\n" 
           . q(    "name": ") . $job . q(",) . "\n" 
           . q(    "modelSeeds": [],) . "\n" 
           . q(    "sequences": [) . "\n" 
           ;

my $footer = q{};
$footer    = q(    ],) . "\n" 
           . q(    "dialect": "alphafoldserver",) . "\n" 
           . q(    "version": 1) . "\n" 
           . q(  }) . "\n"
           . q(]) . "\n"
           ;

$job    = "$job.json";
$job    = safename($job);

foreach my $fasta (@inputs) {
    my $seq = q{};
    open my $FASTA, '<', $fasta;
    while ( my $input = <$FASTA> ) {
        chomp $input;
        if ( $input =~ /\A [>] (\S+) /xms ) {
            $seq = $1;
            $data_ref->{'seq'}->{$seq}->{'res'} = q{};
        }
        elsif ( $input =~ /\S/xms ) {
            $input =~ s/\s//g;
            $data_ref->{'seq'}->{$seq}->{'res'} .= $input;
        }
    }
    close $FASTA;
}

open my $JOB, '>', $job;
my @seqs = sort keys %{ $data_ref->{'seq'} };
print $JOB $header;
my $i = @seqs;
my $j = 0;
foreach my $seq (@seqs) {
    $j++;
    my $res = $data_ref->{'seq'}->{$seq}->{'res'};
    print $JOB q(      {), "\n";
    print $JOB q(        "proteinChain": {), "\n";
    print $JOB q(          "sequence": "), $res, q(",), "\n";
    print $JOB q(          "count": 1), "\n";
    print $JOB q(        }), "\n";
    if ( $i > $j ) {
        print $JOB q(      },), "\n";
    }
    else {
        print $JOB q(      }), "\n";
    }
}
print $JOB $footer;
close $JOB;

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

