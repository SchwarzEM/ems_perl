#!/usr/bin/perl

# extract_mw_data.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/1/2008.
# Purpose: get simple useful CDS/peptide/gene/WBid + pred. MW info table.

use strict;
use warnings;

unless ( $#ARGV == 1 ) { 
    die "Format: ./extract_mw_data.pl wormpep wormpep_EMBOSS_pepstats\n";
}

my $wormpep  = $ARGV[0];
my $pepstats = $ARGV[1];

my $cds     = q{};
my $mwe     = q{};

my %cds2wpe = ();
my %cds2wbg = ();
my %cds2cgc = ();
my %cds2mwe = ();

open my $WORMPEP, '<', $wormpep
    or die "Can't open wormpep file $wormpep: $!";

# Sample wormpep headers:
# >3R5.1 CE24758 WBGene00007065 ...
# >4R79.1a CE35820 WBGene00003525 locus:nas-6 ... 

while (my $input = <$WORMPEP>) { 
    chomp $input;
    if ( $input =~ / \A > 
                        ( \S+ )         # $1 => $cds
                        \s+ 
                        ( CE\d+ )       # $2 => $wpep
                        \s+
                        ( WBGene\d+)    # $3 => $wgene
                        \s+ 
                        locus: ( \S+ )  # $4 => $cgc
                        \s+  
                   /xms ) {
        $cds           = $1;
        my $wpep       = $2;
        my $wgene      = $3;
        my $cgc        = $4;
        $cds2wpe{$cds} = $wpep;
        $cds2wbg{$cds} = $wgene;
        $cds2cgc{$cds} = $cgc;
    } 
    elsif ( $input =~ / \A >
                        ( \S+ )                # $1 => $cds
                        \s+
                        ( CE\d+ )              # $2 => $wpep
                        \s+
                        ( WBGene\d+)           # $3 => $wgene
                        \s+
                   /xms ) {
        $cds           = $1;
        my $wpep       = $2;
        my $wgene      = $3;
        $cds2wpe{$cds} = $wpep;
        $cds2wbg{$cds} = $wgene;
    }
    elsif ( $input =~ /\A>/ ) {  
        die "Couldn't parse header lines of wormpep file $wormpep.\n";
    }
}

close $WORMPEP;

open my $PEPSTATS, '<', $pepstats
    or die "Can't open EMBOSS/pepstats output file $pepstats: $!";

while (my $input = <$PEPSTATS>) { 
    chomp $input;

    # Example input line 1:
    # PEPSTATS of 2L52.1 from 1 to 427

    if ( $input =~ / \A 
                     PEPSTATS 
                     \s+ of \s+ 
                     ( \S+ )                      # $1 => $cds
                     \s+ from \s+ 1 \s+ to \s+ \d+ 
                   /xms ) { 
        $cds = $1;
        $mwe = 0;
    }

    # Example input line 2:
    # Molecular weight = 50017.75             Residues = 427

    elsif ( $input =~ / \A 
                     Molecular \s+ weight \s+ \= \s+ 
                     ( \d+\.\d+ )                  # $1 => $mw 
                   /xms ) { 
        $mwe          = $1;
        if ( !exists $cds2wbg{$cds} ) { 
            die "The protein $cds lacks an official WBGene ID.\n";
        }
        $cds2mwe{$cds} = $mwe;
        $cds            = q{};
        $mwe            = q{};        
    }
}

close $PEPSTATS;

foreach my $cds1 (sort keys %cds2mwe) { 
     print $cds1;
     print "\t";
     print $cds2wpe{$cds1};
     print "\t";
     print $cds2wbg{$cds1};
     print "\t";
     print $cds2cgc{$cds1} if ( exists $cds2cgc{$cds1} );
     print "\t";
     print $cds2mwe{$cds1};
     print "\n";
}

