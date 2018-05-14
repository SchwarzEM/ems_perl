#!/usr/bin/perl

# merge_worm_groups.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/15/2007.
# Purpose: take Boolean gene-trait {0,1} *.csv tables for Cistematic and merge them.

use strict;
use warnings;

my $gene        = "";
my $indiv_trait = "";

my %gene2traits = ();
my @traits      = ();
my @values      = ();
my %obs_traits  = ();

foreach my $infile (@ARGV) { 

    # Import and chomp a file's lines as an array.
    open my $INFILE, $infile or die "Can't open $infile: $!";
    my @inlines = <$INFILE>;
    chomp @inlines;

    # Reset indiv. trait. to null.
    $indiv_trait = "";

    # Shift off the first (header) line and digest it.
    my $header = shift @inlines;

    # The following, and other regexes like it in the program, was hard to get
    #     actually debugged until I figured out what regex Perl would actually follow.
    # '( ... ( ... )? )' or '(?: ...)?' work, but alternative formulations 
    #      like '(?: ...){1,}' -- which should work in Perl -- don't!

    if ( $header !~ / \A WBGene \t \" [^\"]+ (?: \" \t \" [^\t]+)* \" \z /xms ) {
        die "Input file $infile has unparsable header: $header\n";
    }
    if ( $header =~ / \A WBGene \t \" ( [^\"]+ ( \" \t \" [^\t]+)* ) \" \z /xms ) { 
        my $traitline = $1;
        @traits = ();
        @traits = split /\"\t\"/, $traitline;
        foreach $indiv_trait (@traits) { 
            if ( $obs_traits{$indiv_trait} ) { 
                die "$infile specifies a redundant trait: $indiv_trait\n";
            }
            else {
                $obs_traits{$indiv_trait} = 1;
            }
        }
    }
    foreach my $inline (@inlines) { 
        if ( $inline !~ / \A WBGene\d+\S+ \t \d (?: \t \d )* \z /xms ) { 
            die "$infile contains unparsable line: $inline\n";
        }
        if ( $inline =~ / \A ( WBGene\d+\S+ ) ( \t \d ( \t \d )* ) \z /xms ) { 
            $gene = $1;
            my $value_line = $2;
            @values  = ();
            $value_line =~ s/\A\t//;
            @values = split /\t/, $value_line;
            foreach my $i (0..$#values) { 
                $gene2traits{$gene}->{$traits[$i]} = $values[$i];
            }
        }
    }
    close $INFILE;
}

# Print an aggregated header line.
print "WBGene";
foreach $indiv_trait (sort keys %obs_traits) { 
    print "\t\"$indiv_trait\"";
}
print "\n";

foreach my $gene (sort keys %gene2traits) { 
    print "$gene";
    foreach $indiv_trait (sort keys %obs_traits) { 
        if ( $gene2traits{$gene}->{$indiv_trait} ) { 
            print "\t";
            print "$gene2traits{$gene}->{$indiv_trait}";
        }
        else { 
            print "\t";
            print "0";
        }
    }
    print "\n";
}

