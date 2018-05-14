#!/usr/bin/env perl

# get_overall_blastn_matches.pl (formerly, tally_ovis_matches) -- Erich Schwarz <emsch@caltech.edu>, 8/27/2011.
# Purpose: given sequences and their BlastN hits in format-6, compute: size of sequence; amt. of seq. w/ hits; fraction of seq. w/ hits.

use strict;
use warnings;

my $data_ref;

my $seqfile   = $ARGV[0];
my $blastfile = $ARGV[1];

my $seq   = q{};
my $start = q{};
my $stop  = q{};

open my $SEQ, '<', $seqfile or die "Can't open sequence file $seqfile: $!";
while (my $input = <$SEQ>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $seq = $1;
    }
    elsif ( $input =~ / [ACGTNactgn ] /xms ) { 
        $input =~ s/[^ACGTNactgn]//g;
        my $seq_len = length($input);
        $data_ref->{'seq'}->{$seq}->{'length'} += $seq_len;
    }
}
close $SEQ or die "Can't close handle to sequence file $seqfile: $!";

open my $BLAST, '<', $blastfile or die "Can't open BlastN file $blastfile: $!";
while (my $input = <$BLAST>) {
    chomp $input;
    # Hco_McM_000007  OAR24   93.94   66      4       0       71309   71374   20143409        20143474        7e-16    100
    if ( $input =~ /\A (\S+) \t (?: [^\t]* \t){5} (\d+) \t (\d+) /xms ) { 
        $seq   = $1;
        $start = $2;
        $stop  = $3;
        if ( exists $data_ref->{'seq'}->{$seq} ) { 
            foreach my $i ($start..$stop) { 
                $data_ref->{'seq'}->{$seq}->{'marked'}->{$i} = 1;
            }
        }
        if (! exists $data_ref->{'seq'}->{$seq} ) {
            die "Can't recognize sequence in: $input\n";
        }
    }
}
close $BLAST or die "Can't close handle to BlastN file $blastfile: $!";

print "Sequence\tLen.\tMatch\t% match\n";

foreach my $sequence1 ( sort keys %{ $data_ref->{'seq'} } ) { 
    my $seq_len =  $data_ref->{'seq'}->{$sequence1}->{'length'};
    my $match_len = 0;
    if ( exists $data_ref->{'seq'}->{$sequence1}->{'marked'} ) { 
        $match_len = keys %{ $data_ref->{'seq'}->{$sequence1}->{'marked'} };
        my $percent_match = ( ($match_len * 100) / $seq_len);
        $percent_match = sprintf("%.3f", $percent_match);
        my $readable_seq_len   = commify($seq_len);
        my $readable_match_len = commify($match_len);
        print "$sequence1\t$readable_seq_len\t$readable_match_len\t$percent_match\n";
    }
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

