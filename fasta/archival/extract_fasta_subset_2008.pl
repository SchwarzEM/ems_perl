#!/usr/bin/env perl

# extract_fasta_subset_2008.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/24/08.
# Purpose: OLDER version -- given file listing seqs. *or single-name argument*, extract from FASTA file; warn of misses.

use strict;
use warnings;

unless ($#ARGV == 1) { 
    die "Format: extract_fasta_subset.pl  [seqs. list OR seq. name]  [large FASTA to extract]\n";
}

my $input_list     = $ARGV[0];
my $input_fasta    = $ARGV[1];
my %input_names    = ();
my $reading_subset = "no";

my $warnings      = $input_fasta 
                    . "." 
                    . $input_list 
                    . ".warnings"
                    ;

if (-e $input_list) { 
    open (my $INPUT_LIST, "$input_list") 
        or die "Can't open $input_list. $!\n";

    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
            my $seq = $1;
            $input_names{$seq} = 1;
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (!-e $input_list) { 
    $input_names{$input_list} = 1;
}

open my $INPUT_FASTA, '<', $input_fasta 
    or die "Can't open $input_fasta. $!\n";

while (my $fasta_line = <$INPUT_FASTA>) {
    if ( ( $fasta_line =~ / \A > ( \S+ ) \s* /xms ) 
         and ( $input_names{$1} ) ) { 
        my $seq = $1;
        print "$fasta_line";
        $reading_subset = "yes";
        delete $input_names{$seq};  # not 'undef'; Perl Cookbook 5.4.
    }
    elsif ( ( $fasta_line =~ / \A > ( \S+ ) \s* /xms ) 
            and (! $input_names{$1} ) ) {
        $reading_subset = "no";
    }
    elsif ( $reading_subset eq "yes" ) {
        print "$fasta_line";
    }
}
close $INPUT_FASTA 
    or die "Can't close filehandle to FASTA file $input_fasta: $!";

# Scan hash for any non-zero values; if any, warn they were missed.

my $missing_count = scalar (keys %input_names);

if ( $missing_count >= 1) { 
    open my $WARNINGS, '>', $warnings 
        or die "Can't open warnings file $warnings. $!";
    foreach my $key (sort keys %input_names) {
        if ( exists $input_names{$key} ) {   # not 'defined'; Perl Cook. 5.2.
            print {$WARNINGS} "$key not found in $input_fasta\n";
        }
    }
    close $WARNINGS 
        or die "Can't close filehandle to warnings file $warnings: $!";
}

