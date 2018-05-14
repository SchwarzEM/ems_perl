#!/usr/bin/env perl

# extract_integr8_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/5/08.
# Purpose: w/ listed seqs. or 1-name argument, extract subset of Intergr8 .dat; warn of misses.

use strict;
use warnings;

unless ($#ARGV == 1) { 
    die "Format: ./extract_integr8_subset.pl",
        "  [seqs. list OR seq. name]",
        "  [large integr8.dat to extract]\n",
        ;
}

my $input_list     = $ARGV[0];
my $input_integr8  = $ARGV[1];
my %input_names    = ();
my $reading_subset = "no";

my $warnings      = $input_list 
                    . ".warnings"
                    ;

if (-e $input_list) { 
    open (my $INPUT_LIST, "$input_list") 
        or die "Can't open $input_list. $!\n";

    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
            my $protein = $1;
            $input_names{$protein} = 'wanted';
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (!-e $input_list) { 
    $input_names{$input_list} = 'wanted';
}

open my $INPUT_INTEGR8, '<', $input_integr8 
    or die "Can't open $input_integr8. $!\n";

while (my $integr8_line = <$INPUT_INTEGR8>) {

    # N.B.: allow protein record to be read more than once.
    # To change this, second condition should be "eq 'wanted'".
    # 'ID' gives the canonical protein name.

    if ( ( $integr8_line =~ /\A ID \s+ (\S+) /xms ) 
         and ( $input_names{$1} ) ) { 
        my $protein = $1;
        $input_names{$protein} = 'found';        
        print "$integr8_line";
        $reading_subset = 'yes';
    }
    elsif ( ( $integr8_line =~ /\A ID \s+ (\S+) /xms ) 
            and (! $input_names{$1} ) ) {
        $reading_subset = 'no';
    }
    elsif ( $reading_subset eq 'yes' ) {
        print "$integr8_line";
    }
}
close $INPUT_INTEGR8 
    or die "Can't close filehandle to Integr8 file $input_integr8: $!";

# scan hash for any non-zero values, warn that these were missed.

open my $WARNINGS, '>', "$warnings" 
    or die "Can't open warnings file $warnings. $!";
foreach my $key (sort keys %input_names) {
    if ( $input_names{$key} eq 'wanted' ) { 
        print {$WARNINGS} "$key not found in $input_integr8\n";
    }
}
close $WARNINGS 
    or die "Can't close filehandle to warnings file $warnings: $!";
