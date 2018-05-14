#!/usr/bin/env perl

# extract_PFAM_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/21/2011.
# Purpose: given a table of acceptable PFAM hits and a PFAM file (e.g., Pfam-A), extract the subset of HMMs which are in the hitlist.

use strict;
use warnings;
use Getopt::Long;

my $pfam    = q{};
my %opts    = ();
my %motifs  = ();
my $text    = q{};

my $reading = 0;
my $saw_end = 1;

my $input_line_no = 0;

GetOptions(
    'table=s'  => \$opts{'hits'},
    'pfam=s'   => \$opts{'pfam'},
    'negative' => \$opts{'neg'}, 
    'help'     => \$opts{'help'},
);

if ( $opts{'help'} or (! $opts{'hits'}) or (! $opts{'pfam'}) ) { 
    die "Format: extract_PFAM_subset.pl --table|-t [PFAM hits] --pfam|-p [PFAM HMMs file]",
        " --negative|-n [select HMMs *not* in hit table] --help|-h\n",
        ;
}

open my $TABLE, '<', $opts{'hits'} or die "Can't open hit table $opts{'hits'}: $!";
while (my $input = <$TABLE>) { 
    chomp $input;
    if ( $input =~ /\A (?: \S+ \s+){3} (PF\d+\.\d+) \s /xms ) { 
        $pfam          = $1;
        $motifs{$pfam} = 1;
    }
}
close $TABLE or die "Can't close filehandle to hit table $opts{'hits'}: $!";

open my $PFAM, '<', $opts{'pfam'} or die "Can't open PFAM file $opts{'pfam'}: $!";
while (my $input = <$PFAM>) { 
    $input_line_no++;

    # Start reading each motif at the HMMER start line.
    if ( $input =~ /\A HMMER /xms ) { 
        # Should have $reading = 0 every single time.
        if ( ($reading != 0) or ($saw_end != 1) ) { 
            die "Failed to correctly parse data on input line $input_line_no: $input\n";
        }
        else { 
            $reading = 1;
            $saw_end = 0;
            $text .= $input;
        }
    }

    # At the "ACC PFAM" line, decide whether the motif will actually get recorded and printed or not.
    elsif ( $input =~ /\A ACC \s+ (PF\d+\.\d+) \s* \z/xms ) { 
        my $pfam = $1;
        # Should have $reading = 1 every single time.
        if ( ($reading != 1) or ($saw_end != 0) ) { 
            die "Failed to correctly parse data on input line $input_line_no: $input\n";
        }
        elsif ( (  $motifs{$pfam} ) and (! $opts{'neg'} ) ) { 
            $text .= $input;
        }
        elsif ( (! $motifs{$pfam} ) and (  $opts{'neg'} ) ) { 
            $text .= $input;
        }
        else { 
            $text    = q{};
            $reading = 0;
        }
    }
    elsif ( $input =~ / \A \/ \/ \s* \z /xms ) { 
        if ($reading) {
            $text .= $input;
            print $text;   
        }
        $saw_end = 1;
        $text    = q{};
        $reading = 0;
    }
    elsif ($reading) { 
        $text .= $input;        
    }
}
close $PFAM or die "Can't close filehandle to PFAM file $opts{'pfam'}: $!";

# Clear out any straggler text.
print $text if ( $text =~ /\S/xms );

