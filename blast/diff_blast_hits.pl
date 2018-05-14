#!/usr/bin/env perl

# Purpose: print all *actual* hits from positive BLAST report (which, optionally, lack hit in negative BLAST report).

use strict;
use warnings;
use Getopt::Long;

my $reading   = 0;
my $query     = q{};

my $neg_blast = q{};
my $pos_blast = q{};

my %reportable   = ();
my @report_lines = ();
my @output_lines = ();

GetOptions ( 'negative:s' => \$neg_blast,
             'positive:s' => \$pos_blast, );

if ( ! -r $pos_blast ) { 
    die "Format: ./diff_blast_hits.pl",
        " --negative|-n [neg. BLAST] --positive|-p [pos. BLAST]\n",
        ;
} 

if ($neg_blast) { 
    open my $NEG, '<', $neg_blast 
        or die "Can't open negative (filter) BLAST report $neg_blast: $!";
    while ( my $input = <$NEG> ) { 
        if ($input =~ /\A Query= \s+ (\S+) /xms ) { 
            $query = $1;
        } 
        if ( $input =~ / \A \s+ [*]+ \s+ No \s+ hits \s+ found \s+ [*]+ /xms ) { 
            $reportable{$query} = 1;
        }
    }
    close $NEG or die "Can't close filehandle to BLAST report $neg_blast: $!";
}

open my $POS, '<', $pos_blast 
    or die "Can't open positive BLAST report $pos_blast: $!";
while ( my $input = <$POS> ) { 
    if ($input =~ / \A [A-Z]* BLAST [A-Z]* /xms ) { 
        if (@report_lines) { 
            push @output_lines, @report_lines;
            @report_lines = ();
        }
        $reading = 1;
    }
    if ($input =~ /\A Query= \s+ (\S+) /xms ) {
        $query = $1;
        if ( ($neg_blast) and (! $reportable{$query} ) ) { 
             @report_lines = ();
             $reading = 0;
        }
    }
    if ( $input =~ / \A \s+ [*]+ \s+ No \s+ hits \s+ found \s+ [*]+ /xms ) {
        @report_lines = ();
        $reading = 0;
    }
    if ($reading) { 
        push @report_lines, $input;
    }
}
close $POS or die "Can't close filehandle to BLAST report $pos_blast: $!";

if (@report_lines) {
    push @output_lines, @report_lines;
    @report_lines = ();
}

foreach my $output (@output_lines) { 
    print $output;
}

