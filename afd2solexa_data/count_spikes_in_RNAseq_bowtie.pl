#!/usr/bin/env perl

# count_spikes_in_RNAseq_bowtie.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/3/2011.
# Purpose: work around utter lack of usable, documented code for getting real read counts for spikes.

use strict;
use warnings;

# Note that as of 8/3/2011, only four of the following spikes were actually in use!
# 
#     AP2  Sequence Apetala 2 tailed spike
#     Lambda232  Sequence Lambda 9786 tailed spike
#     Lambda11  Sequence Lambda 11300 tailed spike
#     VATG3  Sequence VATG tailed spike

my %spike2size = ( Lambda232 => 9786,
                   Lambda11  => 11934,
                   OBF5      => 1429,
                   EPR1      => 1451,
                   AGP23     => 325,
                   VATG3     => 376,
                   AP2       => 1405,
                   PDF       => 282,   );

# This is apparent order of most to least abundant spikes.
# Expand this list if more spikes are ever added.
my @spike_names = qw ( VATG3
                       Lambda11
                       Lambda232
                       AP2 );

my $mapped_seq = q{};

foreach my $input_file (@ARGV) { 
    my $raw_readcount = 0;
    my %spike2count = ( Lambda232 => 0,
                        Lambda11  => 0,
                        OBF5      => 0,
                        EPR1      => 0,
                        AGP23     => 0,
                        VATG3     => 0,
                        AP2       => 0,
                        PDF       => 0, );

    open my $INPUT_FILE, '<', $input_file or die "Can't open input file $input_file: $!";
    while (my $input = <$INPUT_FILE>) { 
        $raw_readcount++;
        if ( $input =~ /\A [^\t]+ \t [^\t]+ \t (\S+) \t /xms ) { 
            $mapped_seq = $1;
            if ( exists $spike2size{$mapped_seq} ) { 
                $spike2count{$mapped_seq}++;
            }
        }
    } 
    close $INPUT_FILE or die "Can't close filehandle to input file $input_file: $!";

    my $readcount = commify($raw_readcount);
    print "\nTotal reads from $input_file:  $readcount\n\n";

    print "Spike\tReads\tRPKM\n";

    foreach my $spike_seq (@spike_names) { 
        my $reads = commify($spike2count{$spike_seq});
        my $rpkm = ($spike2count{$spike_seq} / ( $raw_readcount * $spike2size{$spike_seq} ) );
        print "$spike_seq\t$reads\t$rpkm\n";
    }
    print "\n";
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

