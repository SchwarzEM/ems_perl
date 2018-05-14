#!/usr/bin/env perl

# filter_func_refinement_01jun2015.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/3/2010.
# Purpose: select only parts of refined FUNC output which have an actual result with a defined p-value threshold.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @input_files = ();

my $overrep;
my $underrep;
my $pvalue = 1;
my $help;

my $header = "GO_term\tp-value\n";

my $data_ref;

GetOptions ( 'input=s{,}' => \@input_files,
             'overrep'    => \$overrep,
             'underrep'   => \$underrep,
             'pvalue=f'   => \$pvalue,
             'help'       => \$help, );

if ( $help or (! @input_files ) or ( $overrep and $underrep ) or ( (! $overrep) and (! $underrep) ) ) { 
    die "Format: filter_func_refinement_01jun2015.pl\n",
        "    --input        [1+ files, or '-' for input stream]\n",
        "    --overrep|-o   [overrepresented GO terms only]\n",
        "    --underrep|-u  [underrepresented GO terms only]\n",
        "    --pvalue|-p    [p-value threshold; default 1]\n",
        "    --help|-h      [print this message]\n",
        ;
}

if ( (! looks_like_number($pvalue) ) or ( $pvalue < 0 ) or ( $pvalue > 1 ) ) {
    die "p-value threshold must be valid number between 0 and 1, not this; $pvalue\n";
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
foreach my $infile (@input_files) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }

    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ / \A 
                         [^\t]+       # GO ontology root description (e.g., "molecular_function")
                         \t ([^\t]+)  # GO term description
                         \t (GO:\d+)  # GO term ID number
                         \t [^\t]*    # + or - ("sign?")
                         \t [^\t]*    # raw_p_underrepresentation_of_variable
                         \t [^\t]*    # raw_p_overrepresentation_of_variable
                         \t ([^\t]+)  # p_underrepresentation_after_refinement
                         \t ([^\t]+)  # p_overrepresentation_after_refinement
                         \t 
                         \z /xms ) {

            my $go_term_desc  = $1;
            my $go_term_id    = $2;
            my $underrep_prob = $3;
            my $overrep_prob  = $4;
            if ($overrep) {
                if ( ( $overrep_prob >= 0 )  and ( $overrep_prob  <= $pvalue ) ) {
                    my $go_term = "$go_term_desc [$go_term_id]";
                    my $go_prob = $overrep_prob;
                    &store_go_term($go_term, $go_prob);
                }
            }
            elsif ($underrep) {
                if ( ( $underrep_prob >= 0 ) and ( $underrep_prob <= $pvalue ) ) {
                    my $go_term = "$go_term_desc [$go_term_id]";
                    my $go_prob = $underrep_prob;
                    &store_go_term($go_term, $go_prob);
                }
            }
        }
        elsif (  ( $input !~ /\A 
                            root_node_name
                            \t node_name
                            \t node_id 
                            \t sign\?
                            \t raw_p_underrepresentation_of_variable = \S+ 
                            \t raw_p_overrepresentation_of_variable = \S+
                            \t p_underrepresentation_after_refinement 
                            \t p_overrepresentation_after_refinement 
                            \z /xms 
                 ) 
              and ( $input !~ / \A 
                            [^\t]+
                            \t [^\t]+
                            \t GO:\d+
                            \t \-
                            \t [^\t]*
                            \t [^\t]*
                            \t [^\t]+
                            \t [^\t]+
                            \t
                            \z /xms
                 ) 
              ) {
            die "Can't parse input line: $input!\n";
        }
    }
    close $INPUT_FILE;
}

my @output_go_terms = sort { $data_ref->{'go_term'}->{$a}->{'go_prob'} 
                             <=> 
                             $data_ref->{'go_term'}->{$b}->{'go_prob'} 
                      }
                      keys %{ $data_ref->{'go_term'} };
foreach my $output_go_term (@output_go_terms) {
    my $output_pvalue = $data_ref->{'go_term'}->{$output_go_term}->{'go_prob'};
    print $header if $header;
    $header = q{};
    print "$output_go_term\t$output_pvalue\n";
}

sub store_go_term {
    my $_go_term = $_[0];
    my $_go_prob = $_[1];
    if ( exists $data_ref->{'go_term'}->{$_go_term} ) {
        die "Redundant record of GO term \"$_go_term\" with p-value $_go_prob\n";
    }
    $data_ref->{'go_term'}->{$_go_term}->{'go_prob'} = $_go_prob;
    return;
}

