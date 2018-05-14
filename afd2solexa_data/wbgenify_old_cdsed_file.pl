#!/usr/bin/env perl

# wbgenify_old_cdsed_file.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/9/2008.
# Purpose: assign WS170/WS186 protein-coding gene IDs to an file w/ CDS gene names (e.g., Ali's).

use strict;
use warnings;
use File::Basename;

my %orig_cds2wbgene = ();
my %wbgene2orig_cds = ();
my %cds2wbgene      = ();
my %wbgene2cds      = ();
my %wbgene2locus    = ();

unless ($#ARGV == 2) {
    die "Format: ./wbgenify_cdsed_file.pl",
        "  wormpep_archival  wormpep_recent",
        "  [cds_named_file.txt|STDIN]\n",
        ;
}

my $wormpep_archival = shift @ARGV;
my $wormpep_recent   = shift @ARGV;

open (my $WORMPEP_ARCHIVAL, "<", "$wormpep_archival") 
    or die "Can't open $wormpep_archival: $!";
    $wormpep_archival = basename($wormpep_archival);
while (my $input = <$WORMPEP_ARCHIVAL>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) \s+ (.+) \z /xms ) {
        my $cds = $1;
        my $header = $2;
        $cds =~ s/[a-z]+\z//;
        if ($header !~ / WBGene\d+ /xms ) {
            die "Couldn't parse in $wormpep_archival: $input\n";
        }
        if ($header =~ / (WBGene\d+) /xms ) {
            my $wbgene = $1;
            $orig_cds2wbgene{$cds}    = $wbgene;
            $wbgene2orig_cds{$wbgene} = $cds;
        }
    }
}
close $WORMPEP_ARCHIVAL;

open (my $WORMPEP_RECENT, "<", "$wormpep_recent")
    or die "Can't open $wormpep_recent: $!";
while (my $input = <$WORMPEP_RECENT>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) \s+ (.+) \z /xms ) {
        my $cds    = $1;
        my $header = $2;
        my $wbgene = q{};
        $cds =~ s/[a-z]+\z//;
        if ($header !~ / (WBGene\d+) /xms ) {
            die "Couldn't parse in $wormpep_recent: $input\n";
        }
        if ($header =~ / (WBGene\d+) /xms ) { 
            $wbgene = $1;
            $cds2wbgene{$cds} = $wbgene;
            $wbgene2cds{$wbgene} = $cds;
        }
        if ($header =~ / locus:(\S+) /xms ) {
            my $locus = $1;
            $wbgene2locus{$wbgene} = $locus;
        }
    }
}
close $WORMPEP_RECENT;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) (\s.*) /xms ) { 
        my $orig_cds  = $1;
        my $remainder = $2;
        my $wbgene = q{};
        if (! $orig_cds2wbgene{$orig_cds} ) { 
            warn "Failed to recognize $orig_cds\n";
            print "WBGeneUnknown|$orig_cds";
        }
        if ( $orig_cds2wbgene{$orig_cds} ) {
            $wbgene = $orig_cds2wbgene{$orig_cds};
            print "$wbgene";
            if ( $wbgene2cds{$wbgene} ) { 
                print "|$wbgene2cds{$wbgene}";
            }
            if (     (! $wbgene2cds{$wbgene}     )
                 and ( $wbgene2orig_cds{$wbgene} ) )  {
                print "|$wbgene2orig_cds{$wbgene}";
                print ":$wormpep_archival";
            }
            if ( $wbgene2locus{$wbgene} ) {
                print "|$wbgene2locus{$wbgene}";
            }
        }
        print "$remainder\n";
    }
}


