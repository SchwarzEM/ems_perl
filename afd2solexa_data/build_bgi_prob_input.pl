#!/usr/bin/env perl

# build_bgi_prob_input.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/28/2011.
# Purpose: given two read/length tables, make a single table which can be fed to bgi_prob.pl.

use strict;
use warnings;
use Getopt::Long;

# Standard input from a given file:

# 	all.reads       uniq.reads      splice.reads    non_multi.reads multi.reads     gene_length
# WBGene00000001|Y110A7A.10|aap-1	122     122     0       122     0       1.500

my $name_col    = q{};
my $reads_col   = q{};
my $len_col     = q{};
my $a_data_file = q{};
my $b_data_file = q{};
my @values      = ();
my $help;

my $data_ref;

GetOptions ( 'a_data_file=s'   => \$a_data_file,
             'b_data_file=s'   => \$b_data_file,
             'name_column:i'   => \$name_col,
             'reads_column:i'  => \$reads_col,
             'length_column:i' => \$len_col,
             'help'            => \$help,           );

# Have human-readable default values:
$name_col  ||= 1;
$reads_col ||= 2;
$len_col   ||= 7;

if ( $help or (! $name_col ) or (! $reads_col ) or (! $len_col ) or (! $a_data_file ) or (! $b_data_file ) ) { 
    die "Format: build_bgi_prob_input.pl\n",
        "           --a_data_file  |-a  [data set A]\n",
        "           --b_data_file  |-b  [data set B]\n",
        "           --name_column  |-n  [default == 1]\n",
        "           --reads_column |-r  [default == 2]\n",
        "           --length_column|-l  [default == 7]\n",
        "\n",
        ;
}

# Reset column values to zero-based Unix:
$name_col--;
$reads_col--;
$len_col--;

my @sources = qw(file_A file_B);

$data_ref->{'file_A'}->{'name'}   = $a_data_file;
$data_ref->{'file_A'}->{'letter'} = 'A';

$data_ref->{'file_B'}->{'name'}   = $b_data_file;
$data_ref->{'file_B'}->{'letter'} = 'B';

foreach my $source (@sources) { 
    open my $DATA, '<', $data_ref->{$source}->{'name'} or die "Can't open data file $data_ref->{$source}->{'letter'} $data_ref->{$source}->{'name'}: $!";
    while (my $input = <$DATA>) {
        chomp $input;
        if ( $input !~ /\A [^\t]+ (?: \t [^\t]*)* \z /xms ) {
            warn "Can't parse input line from data file $data_ref->{$source}->{'letter'} $data_ref->{$source}->{'name'}: $input\n";
        }
        if ( $input =~ /\A [^\t]+ (?: \t [^\t]*)* \z /xms ) { 
            @values = split "\t", $input;
            if (    ( $values[$name_col]  !~ / \S /xms               ) 
                 or ( $values[$reads_col] !~ / \A \d+ \z /xms        ) 
                 or ( $values[$len_col]   !~ / \A \d+ \. \d+ \z /xms ) ) { 
                die "Can't get usable values from $data_ref->{$source}->{'name'} line: $input\n";
            }

            # Enforce one read count per gene:
            if ( exists $data_ref->{'gene'}->{$values[$name_col]}->{'A_readcount'} ) { 
                if ( $values[$reads_col] != $data_ref->{'gene'}->{$values[$name_col]}->{'A_readcount'} ) { 
                    die "Two different readcounts",
                        " ($values[$reads_col] and $data_ref->{'gene'}->{$values[$name_col]}->{'A_readcount'})",
                        " for gene $values[$name_col] in",
                        " data file $data_ref->{$source}->{'letter'} $data_ref->{$source}->{'name'}: $input\n",
                        ;
                }
            }
            $data_ref->{'gene'}->{$values[$name_col]}->{$source}->{'readcount'}  = $values[$reads_col];

            # Keep a running total of all reads seen for the data set:
            $data_ref->{$source}->{'total_reads'}                               += $values[$reads_col];


            # Expand gene sizes from kb to nt, because bgi_prob.pl needs them in nt.
            $values[$len_col] = (1000 * $values[$len_col]);

            # Data should be well enough filtered that most gene lengths are consistent.
            # At this point, if they still aren't, go for longest length observed:
            if ( exists $data_ref->{'gene'}->{$values[$name_col]}->{'length'} ) { 
                if ( $values[$len_col] > $data_ref->{'gene'}->{$values[$name_col]}->{'length'} ) { 
                    $data_ref->{'gene'}->{$values[$name_col]}->{'length'} = $values[$len_col];    
                }
            }
            # If no previous value recorded, record one:
            if (! exists $data_ref->{'gene'}->{$values[$name_col]}->{'length'} ) {
                $data_ref->{'gene'}->{$values[$name_col]}->{'length'} = $values[$len_col];
            }
        }
    }
    close $DATA or die "Can't close filehandle to data file $data_ref->{$source}->{'letter'} $data_ref->{$source}->{'name'}: $!";
}

foreach my $gene (sort keys %{ $data_ref->{'gene'} } ) { 
    print "$gene";
    print "\t";
    print "$data_ref->{'gene'}->{$gene}->{'length'}";
    print "\t";
    print "$data_ref->{'gene'}->{$gene}->{'file_A'}->{'readcount'}";
    print "\t";
    print "$data_ref->{'gene'}->{$gene}->{'file_B'}->{'readcount'}";
    print "\t";
    print "$data_ref->{'file_A'}->{'total_reads'}";
    print "\t";
    print "$data_ref->{'file_B'}->{'total_reads'}";
    print "\n";
}

