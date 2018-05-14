#!/usr/bin/env perl

# extract_fasta_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/11/2011.
# Purpose: given seqlist file or single-name argument, either extract or exclude from FASTA file; optionally, allow matches rather than exact identities; warn of misses.
# N.B. matches are useful for using one list of gene IDs to filter several different sequence files whose members have names that include those gene IDs.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $input_list     = q{};
my @input_fastas   = ();
my $seq            = q{};
my %input_names    = ();  
my %matched_names  = ();
my $reading_subset = 'no';

my $partial;
my $exclude;
my $warnings;
my $help;

my $query         = q{};
my $positives_ref = ();

GetOptions ( 'list=s'     => \$input_list,
             'fasta=s{,}' => \@input_fastas,
             'partial'    => \$partial, 
             'exclude'    => \$exclude,
             'warnings'   => \$warnings,
             'help'       => \$help,     );

if ( (! $input_list) or (! @input_fastas) or $help ) { 
    die "Format: extract_fasta_subset.pl\n",
        "    --list|-l      [seqs. list OR seq. name]\n",
        "    --fasta|-f     [one or more large FASTAs to extract]\n",
        "    --partial|-p   [optional -- allow matches of list names to *parts* of sequence names]\n",
        "    --exclude|-e   [optional -- reject list, keep rest]\n",
        "    --warnings|-w  [optional -- in either positive or negative selection, warn about listed names that were missed]\n",
        "    --help|-h      [print this message]\n",
        ;
}

my $warning_file = basename($input_fastas[0]) . q{.} . basename($input_list) . '.warnings';
$warning_file    = safename($warning_file);

if (-e $input_list) { 
    open my $INPUT_LIST, '<', $input_list or die "Can't open input list $input_list: $!";
    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
            $seq = $1;
            $input_names{$seq} = 1;
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (! -e $input_list) { 
    $input_names{$input_list} = 1;
}

foreach my $input_fasta (@input_fastas) { 
    open my $INPUT_FASTA, '<', $input_fasta or die "Can't open input FASTA $input_fasta: $!\n";
    while (my $fasta_line = <$INPUT_FASTA>) { 
        chomp $fasta_line;
        if ( $fasta_line =~ / \A > ( \S+ ) \s* /xms ) { 
            $seq = $1;
            if (   ( find_match($seq,$partial)      and (! $exclude ) )
                or ( (! find_match($seq,$partial) ) and $exclude      ) ) {
                print "$fasta_line\n";
                $reading_subset = 'yes';
            }
            elsif (   ( find_match($seq,$partial)      and $exclude      ) 
                    or ( (! find_match($seq,$partial) ) and (! $exclude ) ) ) { 
                $reading_subset = 'no';
            }
        }
        elsif ( $fasta_line =~ / \A > /xms ) { 
            die "Can't parse FASTA line: $fasta_line\n";
        }
        # This can only occur to non-header lines after a positive header:
        elsif ( $reading_subset eq 'yes' ) {
            print "$fasta_line\n";
        }
    }
    close $INPUT_FASTA or die "Can't close filehandle to FASTA file $input_fasta: $!";
}

# Optional (because this can pollute output streams) -- for selection of names, warn about any names that were missed:
if ( $warnings ) { 
    my %missed_names = ();
    foreach my $seqname ( sort keys %input_names ) { 
        if (! exists $matched_names{$seqname}) { 
            $missed_names{$seqname} = 1;
        }
    }
    if ( %missed_names ) {
        open my $WARNING_FILE, '>', $warning_file or die "Can't open warning file $warning_file: $!";
        foreach my $key (sort keys %missed_names) {
            print {$WARNING_FILE} "$key was not found in: @input_fastas\n";
        }
        close $WARNING_FILE or die "Can't close filehandle to warning file $warning_file: $!";
    }
}

sub find_match {
    my $_seq     = $_[0];
    my $_partial = $_[1];
    # Look for partial matches of names in list to the given sequence name:
    if ($_partial) {
        foreach my $_seqname ( sort keys %input_names ) {
            # Omit '/xms' to make pattern matches more reliable.
            if ( $_seq =~ /$_seqname/ ) {
                $matched_names{$_seqname} = 1;
                return 1;
            }
        }
        return 0;
    }
    # Or, look for exact matches:
    else { 
        if ( $input_names{$_seq} )  { 
            $matched_names{$_seq} = 1;
            return 1;
        }
        else { 
            return 0;
        }
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


