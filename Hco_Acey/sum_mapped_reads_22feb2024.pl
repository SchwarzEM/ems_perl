#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile     = q{};
my @replicates = ();
my $rep_count  = 0;
my $genes      = 0;
my $data_ref;

my $header = "Genes\tRep_no\tReplicate\tTotal_reads";

$infile = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: sum_mapped_reads_22feb2024.pl [mapped reads table] > [summary of mapped reads per replicate]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;

    my @inputs = split '\t', $input;
    @replicates = @inputs;
    shift @replicates;

    if ( $inputs[0] eq 'Gene' ) {
        # No 'my' before $rep_count here, because we are globally filling in this value from the first row of data columns.
        $rep_count = @replicates;
        $rep_count--;
        foreach my $i (0..$rep_count) {
            $data_ref->{$i}->{'rep'} = $replicates[$i];
        }
    }
    elsif ( $inputs[0] ne 'Gene' ) {
        $genes++;
        foreach my $i (0..$rep_count) {
            $data_ref->{$i}->{'reads'} += $replicates[$i];
        }
    }
}
close $INFILE;

$genes = commify($genes);

foreach my $i (0..$rep_count) {
    my $replicate   = q{};
    my $total_reads = q{};
    my $j = $i;
    $j++;

    $replicate   = $data_ref->{$i}->{'rep'};
    $total_reads = $data_ref->{$i}->{'reads'};
    $total_reads = commify($total_reads);

    print "$header\n" if $header;
    $header = q{};

    print "$genes\t$j\t$replicate\t$total_reads\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

