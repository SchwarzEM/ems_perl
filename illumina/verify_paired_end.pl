#!/usr/bin/env perl

# verify_paired_end.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/28/2011.
# Purpose:  Take four lines at a time (for FASTA input) or eight lines (for default FASTQ input).  Check for superficially correct headers and enforce paired-end-ness.  Check for sequence lines.

use strict;
use warnings;
use Getopt::Long;

my $input = q{};

my %opts  = ();

my $data_ref;
my $read_name = q{};
my $read_stem = q{};

# Things like $opts{'fasta'} make horrible Booleans -- 
#     endless Perl-misery over whether it exists, is defined as 0 or 1, etc.
#     So go for good old simple predefined simple scalars.

my $fasta_input = 0;
my $fastq_input = 0;

# Default suffixes for paired-end read 1 and paired-end read 2:

$opts{'r1'} = '#0\/1';
$opts{'r2'} = '#0\/2';

GetOptions(
    'fasta'             => \$fasta_input,
    'fastq'             => \$fastq_input,
    "r1=s"              => \$opts{'r1'},
    "r2=i"              => \$opts{'r2'}, 
    "help"              => \$opts{'help'},
);

if ( $opts{'help'} or ( $fasta_input and $fastq_input ) ) {
    print "\n";
    print "usage: verify_paired_end.pl\n";
    print "       --fastq      fastQ input (default assumption)\n";
    print "       --fasta      fastA input (optional; mutually incompatible with fastQ)\n";
    print "       --r1         <suffix1>  suffix marking paired-end read 1, default \"#0/1\".\n";
    print "       --r2         <suffix1>  suffix marking paired-end read 2, default \"#0/1\".\n";
    print "       -h|--help    help - print this message\n";
    print "\n";
    exit;
}

if ( (! $fasta_input) and (! $fastq_input) ) { 
    $fastq_input = 1;
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

while (my $input = <>) {
    chomp $input;
    $i++;
    if ( $fastq_input ) {
        $j = ( $i % 8);
    }    
    if ( $fasta_input ) {
        $j = ( $i % 4);
    }    
    if ( $j == 1 ) { 
        # Avoid normal 'xms' parsing to prevent ambiguous regexes.
        if ( $input !~ /\A$lead_pattern\S+/ ) { 
            die "Can't parse putative header: $input\n";
        }

        $read_name = q{};
        if ( $input =~ /\A$lead_pattern(\S+)/ ) {
            $read_name = $1;
        }

        if ( $read_name !~ /\A\S+$opts{'r1'}\z/ ) { 
            die "Can't parse read name $read_name in header: $input\n";
        }
        elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
            die "Header $input somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
        }
        elsif ( $read_name =~ /\A(\S+)$opts{'r1'}\z/ ) {
            # Get the stem name for an r1 read.
            $read_stem = $1;

            # Halt loudly if an identically named r1 has already been read.
            if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                die "Multiple sequence entries with sequence name $read_stem$opts{'r1'}!\n";
            }
            # Halt loudly if its r2 already has been read.
            if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                die "Reads for $read_stem are in wrong order!\n";
            }

            # Record the reading of first (and only allowed) r1 entry.
            $data_ref->{'observed_name'}->{$read_stem}->{'r1'} = 1;
        }
        else { 
            die "Can't parse input: $input\n";
        }
    }

    # Check the FASTA-style sequence associated with namestem/r1 read.
    elsif ( $j == 2 ) { 
        $input =~ s/\s+\z//;
        if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
            die "Can't parse sequence of $read_name: $input\n";
        }
    }


    # On line 5 of 8 (for FASTQ) or line 3 of 4 (for FASTA), require an r2 match for the namestem.
    elsif (    ( $fastq_input and ( $j == 5 ) ) 
            or ( $fasta_input and ( $j == 3 ) ) ) { 

        my $curr_read_stem = q{};
        if ( $input !~ /\A$lead_pattern\S+/ ) {
            die "Can't parse putative header: $input\n";
        }
    
        $read_name = q{};
        if ( $input =~ /\A$lead_pattern(\S+)/ ) {
            $read_name = $1;
        }

        if ( $read_name !~ /\A\S+$opts{'r2'}\z/ ) {
            die "Can't parse read name $read_name in header: $input\n";
        } 
        elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
            die "Header $input somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
        }
        elsif ( $read_name =~ /\A(\S+)$opts{'r2'}\z/ ) {
            # Get the stem name for an r1 read.
            $curr_read_stem = $1;

            # Enforce identical read stems between pairs:
            if ( $curr_read_stem ne $read_stem ) { 
                die "Read stems $read_stem and $curr_read_stem do not match!\n";
            }

            # Halt loudly if an identically named r2 has already been read.
            if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                die "Multiple sequence entries with sequence name $read_stem$opts{'r2'}\n";
            }

            # Halt loudly if its r1 has not been read.
            if (! exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                die "Missing first read for $read_stem!\n";
            }

            # Record the reading of first (and only allowed) r2 entry.
            $data_ref->{'observed_name'}->{$read_stem}->{'r2'} = 1;
        }
        else {
            die "Can't parse input: $input\n";
        }
    }

    # Check the FASTA sequence associated with namestem/r2 read.
    elsif (    ( $fastq_input and ( $j == 6 ) )
            or ( $fasta_input and ( $j == 4 ) ) ) {
        $input =~ s/\s+\z//;
        if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
            die "Can't parse sequence of $read_name: $input\n";
        }
        $read_name = q{};
        $read_stem = q{};
    }
}

#      Use double-quotes so $ARGV actually gets interpolated!

print "The input from $ARGV appears to be a completely correct paired-end ";
print 'FASTQ' if $fastq_input;
print 'FASTA' if $fasta_input;
print " sequence.\n";

