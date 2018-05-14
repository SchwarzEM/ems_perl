#!/usr/bin/env perl

# tally_seqlens.pl -- Erich Schwarz <emsch@caltech.edu>, 8/27/2011.
# Purpose: given a set of sequences in one or more files, tabulate their lengths and note any mismatches.

use strict;
use warnings;

my $seq    = q{};
my $infile = q{};
my $nt_len = q{};

# Allow user to specify order of files:
my @infiles   = @ARGV;
my @sequences = ();

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms) { 
        $seq = $1;
        $infile = $ARGV;
        if ( exists $data_ref->{'infile'}->{$infile}->{'seq'}->{$seq} ) { 
            die "Redundant sequence $seq in input file $infile!\n";
        }
        $data_ref->{'infile'}->{$infile}->{'seq'}->{$seq} = 0;
        $data_ref->{'sequence'}->{$seq} = 1;
    }
    elsif ( $input =~ / [ACGTNacgtn] /xms ) { 
        $input =~ s/[^ACGTNacgtn]//g;
        $nt_len = length($input);
        $data_ref->{'infile'}->{$infile}->{'seq'}->{$seq} += $nt_len;
    }
}

my $header = join "\t", @infiles;
$header = "Sequence\t$header\tLength difference?";
print "$header\n";

@sequences = sort keys %{ $data_ref->{'sequence'} };
foreach my $sequence1 (@sequences) { 
    my $output      = $sequence1;
    my $prev_len    = 'none';
    my $discrepancy = q{};
    foreach my $infile1 (@infiles) { 
        my $length       = 0;
        my $readable_len = 0;
        if ( exists $data_ref->{'infile'}->{$infile1}->{'seq'}->{$sequence1} ) { 
            $length = $data_ref->{'infile'}->{$infile1}->{'seq'}->{$sequence1};
            $readable_len = commify($length);
        }
        $output = $output . "\t$readable_len";
        if ( ( $prev_len ne 'none' ) and ( $prev_len != $length ) ) { 
            $discrepancy = 'Discrepancy';
        }
        $prev_len = $length;
    }
    $output = $output . "\t$discrepancy";
    print "$output\n";
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

