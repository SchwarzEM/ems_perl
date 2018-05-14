#!/usr/bin/perl

# quik_pubmed_reformat.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/14/2007.
# Purpose: quick-and-dirty, but time-saving, reformat of PubMed plaintext ref.

use strict;
use warnings;

my $full_input = q{};

while (my $input = <>) { 
    $input =~ s/\n/ /xms;
    $full_input .= $input;
}
$full_input =~ s/[ ]{2,}/ /g;
$full_input =~ s/\A\s+//;
$full_input =~ s/\s+\z//;

while ( $full_input =~ / \A
                         (.+?) 
                         (\w+) 
                         \s 
                         ([A-Z]+) 
                         , 
                         (.+) 
                         \z
                       /xms ) { 
    my $prequel  = $1;
    my $surname  = $2;
    my $initials = $3;
    my $sequel   = $4;

    $surname .= ","; 
    $initials =~ s/([A-Z])/$1./g; 
    $full_input = $prequel 
                  . $surname 
                  . " " 
                  . $initials 
                  . ","
                  . $sequel
                  ;
}

if ( $full_input =~ / \A
                      (.+)
                      :
                      (\w+)
                      \-
                      (\w+) 
                      \.
                      \z
                    /xms ) {
    my $prequel  = $1;
    my $start_pg = $2;
    my $stop_pg  = $3;

    if ( $prequel =~ / \A 
                       (.+) 
                       (\d+)
                       \( \d+ \) 
                       \z 
                    /xms ) { 
        $prequel = $1 . " " . $2;
    }

    my @start_pp = split //, $start_pg;
    my @stop_pp  = @start_pp;
    my @rev_stop = split //, $stop_pg;

    # '@rev_stop' in scalar context
    foreach my $i (1..@rev_stop) { 
        # overwrite digits
        $stop_pp[-$i] = $rev_stop[-$i];
    } 
    $stop_pg = join q{}, @stop_pp;

    $full_input = $prequel
                  . ", "
                  . $start_pg
                  . "-"
                  . $stop_pg
                  . "."
                  ;
}

print "$full_input\n";

