#!/usr/bin/env perl

# filter_nonhit_BLAST.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/26/2008.
# Purpose: filter BLAST text, omitting non-hit reports.

my $stored_lines = q{};
my $storing      = 0;

my $start_string = '\A [A-Z]* BLAST [A-Z]*';
# my $start_string = '\A BLASTN \s+ 2\.2\.19 \s+ \[Nov\-02\-2008\]';
my $kill_string = '\A \s* [\*]+ \s+ No \s+ hits \s+ found \s+ [\*]+';

while (my $input = <>) { 
    if ( $input =~ / $start_string /xms ) {
        print $stored_lines if $stored_lines;
        $stored_lines = q{};
        $storing  = 1;
        $stored_lines .= $input;
    }
    if ( ($storing) and ( $input !~ / $start_string /xms ) ) {
        $stored_lines .= $input;
    }
    if ( $input =~ / $kill_string /xms ) { 
        $storing  = 0;
        $stored_lines = q{};
    }
}

print $stored_lines if $stored_lines;

