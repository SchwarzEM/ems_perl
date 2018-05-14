#!/usr/bin/env perl

# extract_pfam_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/24/2010.
# Purpose: given motiflist file or single-name argument, either extract or exclude from PFAM file; warn of misses.
# Note: as written, this ducks slight complexity of accepting either NAMEs or ACCs, and just goes with ACCs.

use strict;
use warnings;
use Getopt::Long;

my $input_list = q{};
my $input_pfam = q{};
my $exclude;

my $query         = q{};
my $positives_ref = ();

my $reading_subset  = 1;
my @stored_lines    = ();

GetOptions ( 'list=s'  => \$input_list,
             'pfam=s'  => \$input_pfam, 
             'exclude' => \$exclude,    );

unless ( $input_list and $input_pfam ) { 
    die 'Format: extract_fasta_subset.pl',
        ' --list|-l [seqs. list OR seq. name]',
        ' --pfam|-p [large PFAM to extract]',
        ' --exclude|-e [optional -- reject list, keep rest]',
        "\n",
        ;
}

my %input_names = ();

my $warnings    = $input_pfam 
                  . "." 
                  . $input_list 
                  . ".warnings"
                  ;

if (-e $input_list) { 
    open my $INPUT_LIST, '<', $input_list or die "Can't open $input_list. $!\n";

    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
            my $accession = $1;
            $input_names{$accession} = 1;
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (!-e $input_list) { 
    $input_names{$input_list} = 1;
}

open my $INPUT_PFAM, '<', $input_pfam or die "Can't open $input_pfam. $!";

while (my $pfam_line = <$INPUT_PFAM>) {

# Typical input:
# 
# HMMER3/b [3.0b2 | June 2009]
# NAME  1-cysPrx_C
# ACC   PF10417.2
#   [...]
# //

    # Include the top HMMER line of the accession, but don't rely on it to define reading/motif.
    if ( $pfam_line =~ / \A HMMER\d+.* \z /xms ) { 
        push @stored_lines, $pfam_line;
    }

    # At the end of each PFAM text, print anything stored to print ...
    if ( $pfam_line =~ / \A \/\/ \s* \z /xms ) { 
        if (@stored_lines) { 
            push @stored_lines, $pfam_line;
            print @stored_lines;
        }
        # ... then clear stored lines and stop storing.
        @stored_lines = ();
        $reading_subset = 0;
    }

    # If a clear case occurs of an unwanted ACC, also clear stored lines and stop storing.
    if ( ( $pfam_line =~ / \A (?:ACC) \s+ ( \S+ ) \s* /xms ) 
         and (       ( $exclude      and $input_names{$1}      )
                  or ( (! $exclude ) and (! $input_names{$1} ) )
             )
       ) {
        $reading_subset = 0;
        @stored_lines = ();
    } 

    # Keep track of all *wanted* ACCs, for later error report below:
    if ( ( $pfam_line =~ / \A (?:ACC) \s+ ( \S+ ) \s* /xms ) 
         and (    ( $input_names{$1} and (! $exclude)          ) 
               or ( $exclude         and (! $input_names{$1} ) ) 
             ) 
       ) { 
        my $pfam = $1;
        delete $input_names{$pfam};  # not 'undef'; Perl Cookbook 5.4.
    }

    # If neither at the top of a file nor at the end of an accession...
    if ( ( $pfam_line !~ / \A HMMER\d+.* \z /xms ) and ( $pfam_line !~ / \A \/\/ \s* \z /xms ) ) { 

        # First, see if at the start of an accession; if so, start recording.
        if ( $pfam_line =~ / \A (?:NAME) \s+ \S+ \s* /xms ) {
            $reading_subset = 1;
        }

        # If decided to record (here, or in an earlier line), then store lines.
        if ($reading_subset) { 
            push @stored_lines, $pfam_line;
        }
    }
}
close $INPUT_PFAM or die "Can't close filehandle to PFAM file $input_pfam: $!";

# If this is a 'normal' extraction, not an exclusion-filtering:
if (! $exclude ) { 
    # Scan hash for any non-zero values:
    my $missing_count = scalar (keys %input_names);

    # If any, warn they were missed:
    if ( $missing_count >= 1) { 
        open my $WARNINGS, '>', $warnings 
            or die "Can't open warnings file $warnings. $!";
        foreach my $key (sort keys %input_names) {
            if ( exists $input_names{$key} ) {   # not 'defined'; Perl Cook. 5.2.
                print {$WARNINGS} "$key not found in $input_pfam\n";
            }
        }
        close $WARNINGS 
            or die "Can't close filehandle to warnings file $warnings: $!";
    }
}

