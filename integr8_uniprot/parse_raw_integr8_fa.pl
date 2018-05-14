#!/usr/bin/env perl

# parse_raw_integr8_fa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/4/2008.
# Purpose: for huge FASTA @ indiv. Integr8 proteomes, make names nonredundant.

use strict;
use warnings;

my %prot_header     = ();
my %printed         = ();
my $ok_to_print_seq = 0;

if ($#ARGV > 0) { 
    die "Format: ./parse_raw_integr8_fa.pl  [raw_integr8.fa]\n";
}

# First opening of raw file.
my $input_file = $ARGV[0];
open my $RAW_INT8, '<', $input_file 
    or die "Can't open putative raw Integr8 FASTA file $input_file.\n";

# First pass: record names and headers non-redundantly; track multiples.
while (my $input = <$RAW_INT8>) { 
    chomp $input;
    if ($input =~ /\A > (\S+) (.*) \z /xms ) { 
        my ($name, $header);
        $name   = $1;
        $header = $2;
        if ( $prot_header{$name} ) { 
            $header = q{|} . $header;
        }
        $prot_header{$name} .= $header;
    }
}
close $RAW_INT8 
    or die "Failed to close filehandle for input file $input_file: $!";

# Second opening of raw file:
open $RAW_INT8, '<', $input_file 
    or die "Can't open putative raw Integr8 FASTA file $input_file.\n";

# Second pass: this time, print each named sequence once w/ all headers.
while (my $input = <$RAW_INT8>) {
    chomp $input;
    if (     ( $input !~ /\A > (\S+) (.*) \z /xms )
         and ( $ok_to_print_seq                   ) ) { 
        print "$input\n";
    }
    if ( $input =~ /\A > (\S+) /xms ) {
        my $name;
        $name  = $1;
        if ( exists $printed{$name} ) { 
            $ok_to_print_seq = 0;
        }
        if (! exists $printed{$name} ) { 
            $ok_to_print_seq = 1;
            print '>',
                  $name,
                  $prot_header{$name},
                  "\n",
                  ;
            $printed{$name} = 1;
        }
    }
}
close $RAW_INT8
    or die "Failed to close filehandle for input file $input_file: $!";

