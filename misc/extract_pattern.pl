#!/usr/bin/env perl

# extract_pattern.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/14/2009.
# Purpose: extract string- or file-defined pattern from text.

# Example 'pattern' file (must be one line and must have parentheses):
# WB:(\[.+\])\t

use strict;
use warnings;
use Getopt::Long;

my $pattern;
my $pattern_file;
my $input_file;
my $output_file;

my $usage_string = "Syntax: ./extract_pattern.pl" 
                   . " ( -p 'pattern_string' [need '...' for bash]"
                   . " | -f pattern_file )" 
                   . " -i input_file ( STDOUT | -o output_file )"
                   . "\n"
                   ;

GetOptions( "pattern=s" => \$pattern,
            "file=s"    => \$pattern_file,
            "input=s"   => \$input_file,
            "output=s"  => \$output_file, );

if ( ($input_file) and (! -r $input_file ) ) { 
    warn "$usage_string";
    die "Can't read input file $input_file.\n";
}

if ( ( (! $pattern ) or ( $pattern !~ /\S/xms ) ) 
     and ( (! $pattern_file ) or (! -r $pattern_file ) ) ) { 
    warn "$usage_string";
    if (! $pattern =~ /\S/xms ) {
        warn "No pattern string given.\n";
    }
    if ( (! $pattern_file ) or (! -r $pattern_file ) ) {
        warn "No readable pattern file given.\n";
    }
    die "Must provide either \"-p 'pattern_string'\" or",
        " \"-f pattern_file\" argument.\n",
        ;
}

if ( ( ($pattern) and ( $pattern =~ /\S/xms ) ) 
     and ( ($pattern_file) and (-r $pattern_file ) ) ) { 
    warn "Ignoring readable pattern file $pattern_file",
         " in preference for pattern string $pattern!\n",
         ;
}

if  ( ( ($pattern_file) and (-r $pattern_file ) ) and (! $pattern ) ) { 
    open my $PATTERN, '<', $pattern_file 
        or die "Cannot open $pattern_file: $!\n";
    my @pattern_lines = <$PATTERN>;
    chomp @pattern_lines;
    $pattern = $pattern_lines[0];
    close $PATTERN;
}

if ( $pattern !~ / \( \S+ \) /xms ) {
    warn "$usage_string";
    die "Pattern $pattern lacks required string-capturing parentheses.\n";
}

open my $INPUT, '<', $input_file   or die "Can't open $input_file $!\n";

my $OUTPUT;
if ($output_file) { 
    open $OUTPUT, '>', $output_file or die "Can't open $output_file $!\n";
    # Redirect STDOUT to $OUTPUT!
    select $OUTPUT;
}

while (my $input = <$INPUT> ) {
    chomp $input;
    # N.B.: $pattern MUST include '\(.+\)'.
    if ($input =~ /$pattern/) {
        # 'print' goes to $OUTPUT if $output was specified.
        print "$1\n";
    }
}

close $INPUT;
if ($output_file) { 
    close $OUTPUT;
}

