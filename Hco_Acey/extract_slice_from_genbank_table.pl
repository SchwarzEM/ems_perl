#!/usr/bin/env perl

# Erich Schwarz <ems@emstech.org>, 2/21/2014.
# Purpose: given a list of one or more sequence names and a feature table intended for GenBank Sequin, extract only those parts of the table annotating those sequences.

use strict;
use warnings;

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $input_list     = q{};
my @input_tables   = ();
my %input_names    = ();  
my %matched_names  = ();
my $reading_subset = 'no';

my $exclude;
my $warnings;
my $help;

GetOptions ( 'list=s'     => \$input_list,
             'table=s{,}' => \@input_tables,
             'exclude'    => \$exclude,
             'warnings'   => \$warnings,
             'help'       => \$help,     );

if ( (! $input_list) or (! @input_tables) or $help ) { 
    die "Format: extract_slice_from_genbank_table.pl\n",
        "    --list|-l      [seqs. list OR seq. name]\n",
        "    --table|-t     [one or more large NCBI feature tables to extract]\n",
        "    --exclude|-e   [optional -- reject list, keep rest]\n",
        "    --warnings|-w  [optional -- in either positive or negative selection, warn about listed names that were missed]\n",
        "    --help|-h      [print this message]\n",
        ;
}

my $warning_file = basename($input_tables[0]) . q{.} . $input_list . '.warnings';
$warning_file    = safename($warning_file);

if (-e $input_list) { 
    open my $INPUT_LIST, '<', $input_list or die "Can't open input list $input_list: $!";
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
if (! -e $input_list) { 
    $input_names{$input_list} = 1;
}

foreach my $input_table (@input_tables) { 
    open my $INPUT_TABLE, '<', $input_table or die "Can't open input table $input_table: $!\n";
    while (my $table_line = <$INPUT_TABLE>) { 
        chomp $table_line;
        if ( $table_line =~ / \A > Features \s+ (\S+) \b /xms ) { 
            my $seq = $1;
            $reading_subset = 'no';
            if ( ( $input_names{$seq} and (! $exclude ) ) 
                 or 
                 ( (! $input_names{$seq} ) and $exclude ) 
               ) {
                print "$table_line\n";
                $matched_names{$seq} = 1;
                $reading_subset = 'yes';
            }
        }
        # We have to allow an exception for NCBI's shoddy syntax, which can easily have '>\d+' at the start of a line.
        elsif ( ( $table_line =~ / \A > /xms ) and ( $table_line !~ / \A > \d+ /xms ) ) { 
            die "Can't parse table line: $table_line\n";
        }
        # This can only occur to non-header lines after a positive header:
        elsif ( $reading_subset eq 'yes' ) {
            print "$table_line\n";
        }
    }
    close $INPUT_TABLE or die "Can't close filehandle to input table $input_table: $!\n";
}

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
            print {$WARNING_FILE} "$key was not found in: @input_tables\n";
        }
        close $WARNING_FILE or die "Can't close filehandle to warning file $warning_file: $!";
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

