#!/usr/bin/env perl

# count_mispaired_ends.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/10/2012.
# Purpose:  On a file that is *supposed* to be paired-end, take four or eight lines at a time (depending on FASTA or (default) FASTQ input).  Check for superficially correct headers and DNA sequence.  Count mispaired reads.  Defaults to new-fangled '[ ]1' and '[ ]2' read markers, not older '#0\/1' and '#0\/2'.

use strict;
use warnings;
use Getopt::Long;

my $input        = q{};
my $read_name    = q{};
my $read_stem    = q{};
my $jumbled_pair = 0;

my $data_ref;

my $ordered_count   = 0;
my $jumbled_count   = 0;
my $total_count     = 0;
my $percent_jumbled = 0;

my %opts  = ();

my @input_seqs  = ();
my $fasta_input = 0;
my $fastq_input = 0;

# Default suffixes for paired-end read 1 and paired-end read 2:
$opts{'r1'} = '[ ]1';
$opts{'r2'} = '[ ]2';

# Note that these used to be the older-style defaults:
# $opts{'r1'} = '#0\/1';
# $opts{'r2'} = '#0\/2';
# But those don't work well as typed arguments, so enter " --r1 [#]0\/1 --r2 [#]0\/2 " instead.

GetOptions(
    'input=s{,}' => \@input_seqs,
    'fasta'      => \$fasta_input,
    'fastq'      => \$fastq_input,
    "r1=s"       => \$opts{'r1'},
    "r2=s"       => \$opts{'r2'}, 
    "help"       => \$opts{'help'},
);

if (    (! @input_seqs) 
     or (! $opts{'r1'}   ) 
     or (! $opts{'r2'}   ) 
     or ( $fasta_input and $fastq_input ) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) { 
    print "\n";
    print "usage: count_mispaired_ends.pl\n";
    print "       -i|--input            <name1+>   input one or more file names, or \"-\" for STDIN stream; must be all fasta or fastq format; mandatory\n";
    print "       --fastq                          fastQ format for input and outputs (default)\n";
    print "       --fasta                          fastA format for input and outputs (mutually exclusive with fastQ)\n";
    print "       --r1                  <suffix1>  suffix marking paired-end read 1, default \"[ ]1\"; for old-style reads, type \"", '[#]0\/1', "\".\n";
    print "       --r2                  <suffix2>  suffix marking paired-end read 2, default \"[ ]2\"; for old-style reads, type \"", '[#]0\/2',  "\".\n";
    print "       -h|--help                        help - print this message\n";
    print "\n";
    exit;
}

if ( (! $fasta_input ) and (! $fastq_input ) ) { 
    $fastq_input = 1;
}

my $INFILE;

foreach my $input_seq (@input_seqs) { 
    # Accept either a stream from '-' or a standard file.
    if ($input_seq eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE, '<', $input_seq or die "Cannot open input file $input_seq: $!";
    }

    my $i = 0;
    my $j = 0;
    my $lead_pattern = q{};
    if ( $fastq_input ) {
        $lead_pattern = '[@]';
    }
    if ( $fasta_input ) {
        $lead_pattern = '>';
    }

    while (my $input = <$INFILE>) {
        chomp $input;
        $i++;
        $j = correct_modulus($i);

        if ($j == 1) { 
            $jumbled_pair = 0; 
            $read_name = q{};
            $read_stem = q{};

            # Avoid normal 'xms' parsing to prevent ambiguous regexes.
            # Enforce deterministic parsing; disallow silent parse-failure.

            if ( $input !~ /\A$lead_pattern\S+/ ) { 
                die "ERROR 1: Can't parse lead pattern of putative header from first member of pair: $input\n";
            }
            elsif ( $input =~ /\A$lead_pattern(\S+(?:$opts{'r1'}|$opts{'r2'}))/ ) {
                $read_name = $1;
                # Nested round of deterministic parsing:
                if ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
                    die "ERROR 2: Indeterminate read name from first member of pair: $read_name\n";
                }
                elsif ( $read_name =~ /\A(\S+)$opts{'r1'}\z/ ) {
                    # Get the stem name for an r1 read.
                    $read_stem = $1;
                    $data_ref->{'observed_name'}->{$read_stem}->{'r1'} = 1;
                }
                elsif ( $read_name =~ /\A(\S+)$opts{'r2'}\z/ ) {
                    $read_stem = $1;
                    $jumbled_pair = 1;
                }
                else {
                    die "ERROR 3: Can parse lead pattern and --r1/--r2, but otherwise can't parse input from first member of pair: $input\n";
                }
            }
            else { 
                die "ERROR 4: Can parse lead pattern but not --r1/--r2 of input from first member of pair: $input\n";
            }
        }

        # Check the FASTA sequence associated with namestem/r1 read.
        elsif ( $j == 2 ) { 
            $input =~ s/\s+\z//;
            if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
                die "ERROR 5: Can't parse sequence of $read_name from first member of pair: $input\n";
            }
        }

        # On line 5/8 or 3/4, require an r2 match for the namestem.
        elsif (    ( $fastq_input and ($j == 5) ) 
                or ( $fasta_input and ($j == 3) ) ) {
            my $curr_read_stem = q{};
            $read_name = q{};
            if ( $input !~ /\A$lead_pattern\S+/ ) {
                die "ERROR 6: Can't parse lead pattern of putative header from second member of pair: $input\n";
            }
            elsif ( $input =~ /\A$lead_pattern(\S+(?:$opts{'r1'}|$opts{'r2'}))/ ) {
                $read_name = $1;
                if ( $read_name !~ /\A\S+$opts{'r2'}\z/ ) {
                    $jumbled_pair = 1;
                } 
                elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
                    die "ERROR 7: Indeterminate read name from second member of pair: $read_name\n";
                }
                elsif ( $read_name =~ /\A(\S+)$opts{'r2'}\z/ ) {
                    # Get the stem name for an r1 read.
                    $curr_read_stem = $1;

                    # Enforce identical read stems between pairs:
                    if ( $curr_read_stem ne $read_stem ) { 
                        $jumbled_pair = 1;
                    }

                    # Enforce strict order of read pairs:
                    elsif (! exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                        $jumbled_pair = 1;
                    }
                    else { 
                        $jumbled_pair = 0;
                        delete $data_ref->{'observed_name'}->{$read_stem}->{'r1'};
                    }
                }
            }
            else {
                die "ERROR 8: Can parse lead pattern, but can't parse --r2 in input from second member of pair: $input\n";
            }
        }

        # On line 6/8 or 4/4, check the FASTA sequence associated with namestem/r2 read.
        elsif (    ( $fastq_input and ($j == 6) )
                or ( $fasta_input and ($j == 0) ) ) {  # *Not* $j == 4! because 4/4 in Perl modulus = 0!
            $input =~ s/\s+\z//;
            if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
                die "ERROR 9: Can't parse sequence of $read_name from second member of pair: $input\n";
            }
            $read_name = q{};
            $read_stem = q{};
            if ($jumbled_pair) {
                $jumbled_count++;
            } 
            elsif (! $jumbled_pair) {
                $ordered_count++;
            }
        }
    }
    close $INFILE or die "Cannot close filehandle to input file $input_seq: $!";
}

$total_count = $ordered_count + $jumbled_count;

if ($total_count >= 1) {
    $percent_jumbled = 100 * ( $jumbled_count / $total_count );
}
$percent_jumbled = sprintf("%.2f", $percent_jumbled);

$ordered_count = commify($ordered_count);
$jumbled_count = commify($jumbled_count);
$total_count   = commify($total_count);

print "Total pairs:       $total_count\n";
print "Ordered pairs:     $ordered_count\n";
print "Jumbled pairs:     $jumbled_count\n";
print "Percent jumbled:   $percent_jumbled%\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

sub correct_modulus { 
    my $_i = $_[0];
    my $_j = 0;
    if ( $fastq_input ) {
        $_j = ($_i % 8);
    }
    
    if ( $fasta_input ) {
        $_j = ($_i % 4);
    }
    return $_j;
}

