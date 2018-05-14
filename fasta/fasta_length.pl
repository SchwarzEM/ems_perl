#!/usr/bin/env perl

# fasta_length.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/29/2008.
# Purpose: censors too short (e.g., zero-length) or too long FASTAs.

use strict;
use warnings;
use Getopt::Long;

my $min;
my $max;
my $help;

GetOptions ( 'min=i' => \$min,
             'max=i' => \$max, 
             'help'  => \$help, );

if ($help) { 
    die "Format: fasta_length.pl --min [non-negative integer] --max [positive integer] -h|--help\n";
}

# Error-check minimum and maximum (if any were provided):

if ( (defined $min) and ( ( $min < 0 ) or ( $min != int($min) ) ) ) { 
    die "Minimum value $min should be a nonnegative integer!\n";
}

if ( (defined $max) and ( ( $max < 1 ) or ( $max != int($max) ) ) ) {
    die "Maximum value $max should be a positive integer!\n";
}
if ( (defined $min) and (defined $max) and ($min > $max) ) { 
    die "Minimum value $min should be less than or equal to",
        " maximum value $max!\n",
        ;
}

my $output_line     = q{};
my @output_lines    = ();
my $seq_name        = q{};
my %seqs2headers    = ();
my %sequences       = ();

# Record the incoming FASTA data.

while (my $input_line = <>) { 
    chomp $input_line;
    if ($input_line =~ /\A > ( (\S+) .*) /xms) { 
        $seq_name = $2; 
        $seqs2headers{$seq_name} = $1;
        $sequences{$seq_name} = q{};
    }
    elsif ( $input_line =~ /\S/xms ) { 
        $sequences{$seq_name} .= $input_line;
    }
}

# Weed out bad sequences; print out the good ones.

foreach my $seq_name2 (sort keys %sequences) { 
    my $length = length ( $sequences{$seq_name2} );
    if (  (  (! defined $min ) 
             or 
             ( ( defined $min ) and ( $length >= $min ) )  ) 
          and 
          (  (! defined $max )
             or
             ( ( defined $max ) and ( $length <= $max ) )  ) 
       ) {
        print ">$seqs2headers{$seq_name2}\n";
        @output_lines
            = unpack( "a60" 
              x (length($sequences{$seq_name2})/60 + 1), 
              $sequences{$seq_name2} );
        foreach $output_line (@output_lines) {
            if ($output_line =~ /\S/) {
                    print "$output_line\n";
            }
        }
    }
}

