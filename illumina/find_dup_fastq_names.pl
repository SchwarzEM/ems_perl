#!/usr/bin/env perl

# find_dup_fastq_names.pl -- Erich Schwarz <emsch@caltech.edu>, 9/2/2011.
# Purpose: given a stream of FASTQ text, check whether read names recur (a test of redundancy) and give readcounts.

use strict;
use warnings;

my $i            = 0;
my $readname     = q{};
my %first_seen   = ();
my %repeated     = ();
my %file_w_first = ();
my %file_w_rep   = ();

my $list_of_inputs = join '; ', @ARGV;

while (my $input = <>) { 
    chomp $input;
    $i++;
    if (($i % 4) == 1) { 
        # Sample header:
        # @ILLUMINA-33A494_0009:8:1:2750:1211#0/1
        if ( $input =~ /\A [@] (\S+) /xms ) { 
            $readname = $1;
            if ( exists $first_seen{$readname} ) { 
                $repeated{$readname} = 1;
                $file_w_rep{$ARGV}   = 1;
            }
            else { 
                $first_seen{$readname} = 1;
                $file_w_first{$ARGV}   = 1;
            }
        }
        else { 
            die "Can't parse input from $ARGV [line $i of stream]: $input\n";
        }
    }
}

my $read_count = 0;
if ( keys %first_seen ) { 
    $read_count = keys %first_seen;
}
my $repeat_count = 0;
if ( keys %repeated ) { 
    $repeat_count = keys %repeated;
}
$read_count   = commify($read_count);  
$repeat_count = commify($repeat_count);

my @files_with_first_reads = ();
if ( keys %file_w_first ) { 
    @files_with_first_reads = keys %file_w_first;
}
my @files_with_repeat_reads = ();
if ( keys %file_w_rep ) { 
    @files_with_repeat_reads = keys %file_w_rep;
}
my $list_of_first_read_files  = join q{, }, @files_with_first_reads;
my $list_of_repeat_read_files = join q{, }, @files_with_repeat_reads;

print "\n";
print "Read name count:     $read_count\n";
print "Repeated name count: $repeat_count\n";
print "\n";
print "Files with unique reads:   $list_of_first_read_files\n";
print "Files with repeated reads: $list_of_repeat_read_files\n";
print "\n";


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

