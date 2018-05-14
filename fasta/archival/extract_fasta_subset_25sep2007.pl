#!/usr/bin/perl

# extract_fasta_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/25/07.
# Purpose: given listed sequence set, extract it from FASTA file; warn about misses.

use strict;
use warnings;

unless ($#ARGV == 1) { 
    die "Format: extract_fasta_subset.pl  [seqs. list]  [large FASTA to extract]\n";
}

my $input_list     = $ARGV[0];
my $input_fasta    = $ARGV[1];
my %input_names    = ();
my $reading_subset = "no";

my $output_subset = $input_fasta 
                    . "." 
                    . $input_list 
                    . ".subset"
                    ;

my $warnings      = $input_fasta 
                    . "." 
                    . $input_list 
                    . ".warnings"
                    ;

open (my $INPUT_LIST, "$input_list") 
    or die "Can't open $input_list. $!\n";

while (my $inline = <$INPUT_LIST>) { 
    chomp $inline;
    if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
        $input_names{$1} = 1;
    }
}
close $INPUT_LIST;

open (my $INPUT_FASTA, "$input_fasta") 
    or die "Can't open $input_fasta. $!\n";

while (my $fasta_line = <$INPUT_FASTA>)
{
    if ( ( $fasta_line =~ / \A > ( \S+ ) \s* /xms ) 
         and ( $input_names{$1} ) ) {
        print "$fasta_line";
        $reading_subset = "yes";
        $input_names{$1} = 0;
    }
    elsif ( ( $fasta_line =~ / \A > ( \S+ ) \s* /xms ) 
            and (! $input_names{$1} ) ) {
        $reading_subset = "no";
    }
    elsif ( $reading_subset eq "yes" ) {
        print "$fasta_line";
    }
}
close $INPUT_FASTA ;

# scan hash for any non-zero values, warn that these were missed.
open (my $WARNINGS, ">$warnings") or die "Can't open $warnings. $!\n";
foreach my $key (sort keys %input_names) {
    if ($input_names{$key} == 1) { 
        print $WARNINGS "$key not found in $input_fasta\n";
    }
}
close $WARNINGS;
