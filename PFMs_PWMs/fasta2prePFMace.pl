#!/usr/bin/env perl

# fasta2prePFMace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/15/2010.
# Purpose: Given a pure or FASTA seqlist, print usable proto_PFM.ace.

use strict;
use warnings;
use File::Basename;
use TFBS::PatternGen::SimplePFM;

# Initialize variables:
my $name     = 'default_name';
my $filename = 'default_filename';

if ( ($ARGV[0]) and (-e $ARGV[0]) ) { 
    $filename = basename($ARGV[0]);
    $name = $filename . '_protoPFM.ace';
}

my @sequences  = ();
my $site_count = 0;

# Read in sequence lines from raw list or FASTA:
while (my $input = <>) { 
    chomp $input;
    if ($input !~ /\A >/xms) { 
        $input =~ s/\s//g;
        if ( $input =~ /[^ACGTacgt]/xms ) { 
            die "Can't parse: $input\n";
        }
        if ( $input =~ / \A [ACGTacgt]+ \z /xms ) {
            push @sequences, $input;
        }
    }
}

$site_count = @sequences;

# Print header:
print "Position_Matrix : \"$name\"\n";
print "Description       \"Generated by TFBS::PatternGen::SimplePFM from $filename.\"",
                         " Paper_evidence \"WBPaperXXXXXXXX\"  // pmidXXXXXXX\n",
                         ;
print "Type              Frequency\n";

# Generate, format, and print frequency lines:
my $pfm = TFBS::PatternGen::SimplePFM->new(-seq_list=>\@sequences);
my @countlines = split /\n/, ($pfm->pattern->prettyprint);

foreach my $line (@countlines) { 
    $line =~ s/(\[|\])//g;
    $line =~ s/\A[ ]+//;
    $line =~ s/[ ]+\z//;
    $line =~ s/[ ]{2,}/ /g;
    print "Site_values       $line\n";
}

# This will not work with pre-Jan. 2010 .ace motif models:
print "Sites_used       $site_count\n";


