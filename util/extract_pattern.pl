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
my $output_file;
my $simple = 0;
my @pattern_lines = ();
my $help;

my $usage_string = "Syntax: ./extract_pattern.pl" 
                   . " ( -p 'pattern_string' [need '...' for bash]"
                   . " | -f pattern_file )" 
                   . " -s [option to do simple matches]"
                   . "\n"
                   ;

GetOptions( "pattern=s" => \$pattern,
            "file=s"    => \$pattern_file,
            "output=s"  => \$output_file, 
            "simple=i"  => \$simple,
            "help"      => \$help          );

if ( $help 
     or  (     ( (! $pattern      ) or ( $pattern !~ /\S/xms ) ) 
           and ( (! $pattern_file ) or (! -r $pattern_file   ) ) ) 
          ) { 
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

if (     ( ($pattern     ) and ( $pattern =~ /\S/xms ) ) 
     and ( ($pattern_file) and (-r $pattern_file     ) ) ) { 
    warn "Ignoring readable pattern file $pattern_file",
         " in preference for pattern string $pattern!\n",
         ;
}

if  ( ( ($pattern_file) and (-r $pattern_file ) ) and (! $pattern ) ) { 
    open my $PATTERN, '<', $pattern_file 
        or die "Cannot open $pattern_file: $!\n";
    @pattern_lines = <$PATTERN>;
    chomp @pattern_lines;
    close $PATTERN;
}

my $OUTPUT;
if ($output_file) { 
    open $OUTPUT, '>', $output_file or die "Can't open $output_file $!\n";
    # Redirect STDOUT to $OUTPUT!
    select $OUTPUT;
}

while (my $input = <> ) {
    chomp $input;
    # N.B.: unless opting for 'simple, $pattern MUST include '\(.+\)'.
    # This has been rewritten to try to look for more than one line in a pattern file:

    if (! $simple) { 
        if (@pattern_lines) { 
            foreach my $pat1 (@pattern_lines) { 
                reality_check($pat1);
                if ( $input =~ /$pat1/ ) {
                    # 'print' goes to $OUTPUT if $output was specified.
                    print "$1\n";
                }
            }
        }
        if (! @pattern_lines) { 
            if ( $input =~ /$pattern/ ) {
                reality_check($pattern);
                # 'print' goes to $OUTPUT if $output was specified.
                print "$1\n";
            }
        }
    }
    if ( ($simple) and ( $input =~ /$pattern/ ) ) { 
        print "$input\n";
    }
}

if ($output_file) { 
    close $OUTPUT;
}

sub reality_check { 
    my $pattern_line = $_[0];
    if ( $pattern_line !~ / \( .+ \) /xms ) {
        warn "$usage_string";
        die "Pattern $pattern lacks required string-capturing parentheses.\n";
    }
    return;
}

