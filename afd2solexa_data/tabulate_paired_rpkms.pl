#!/usr/bin/env perl

# tabulate_paired_rpkms.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/23/2008.
# Purpose: get an Excel-parsable table from RPKM files.

use strict;
use warnings;
use File::Basename;

# WBGene -> filename -> RPKM value.
my $gene2values_ref;
my %seen;

foreach my $infile (@ARGV) { 
     my $filename = basename($infile);
     if ( $seen{$filename} ) { 
         die "Redundant file name $filename\n";
     }
     $seen{$filename} = 1;
     open my $INFILE, '<', $infile or die "Can't open $infile: $!";
     while (my $input = <$INFILE>) { 
         chomp $input;
         if ( $input =~ /\A (\S+) \s+ \S+ \s+ (\S+) /xms ) { 
             my $cds  = $1;
             my $rpkm = $2;
             $gene2values_ref->{$cds}->{$filename} = $rpkm;
         }
     }
     close $INFILE or die "Can't close filehandle to $infile: $!";
}

# Print table header:
print "Gene";
foreach my $source1 (sort keys %seen) { 
    print "\t$source1";
}
print "\n";

# Print table:
foreach my $cds_seen (sort keys %{ $gene2values_ref } ) { 
    print $cds_seen;
    foreach my $source2 (sort keys %seen) {
        if ( exists $gene2values_ref->{$cds_seen}->{$source2} ) { 
            print "\t";
            print $gene2values_ref->{$cds_seen}->{$source2};
        }
        else { 
            print "\t0.00";
        }
    }
    print "\n";
}

