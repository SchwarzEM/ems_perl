#!/usr/bin/env perl

# sieve_bed.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/7/2009.
# Purpose: given FASTA winfo headers, BlastN bad hits, Ali's BED -> print bad- or non-bad BED lines.

use strict;
use warnings;
use Getopt::Long;

my %seen_elements = ();
my %blastn_sites  = ();

my $badprint   = 1;
my $goodprint;
my $fasta      = q{};
my $blastn     = q{};
my $bed        = q{};

GetOptions ( 'fasta=s'  => \$fasta,
             'nblast=s' => \$blastn, 
             'bed=s'    => \$bed, 
             'misses'   => \$goodprint, );

if ( (! $fasta) or (! $blastn) or (! $bed) or ( (! $badprint) and (! $goodprint) ) ) {
    die "Format: sieve_bed.pl --fasta|-fa [FASTA file] --nblast|-n ['bad' BlastN file ]",
        " --bed|-b [BED file] --misses|-m [boolean: print the misses]\n",
        ;
}

if ($goodprint) { 
    $badprint = 0;
}

open my $FASTA, '<', $fasta or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) { 
    my $name     = q{};
    my $site     = q{};
    my $chr      = q{};
    my $start_nt = q{};
    my $stop_nt  = q{};

# Sample input:
# >ce6_ct_phastv3filtered_chrIV.8293 range=chrIV:4930944-4931058 5'pad=0 3'pad=0 strand=+ repeatMasking=none

    if ( $input =~ /\A > (\S+) .+ range= (chr \w+) : (\d+) \- (\d+) /xms ) { 
        $name     = $1;
        $chr      = $2;
        $start_nt = $3;
        $stop_nt  = $4;
        $start_nt-- ;    # Map from FASTA numbering to BED numbering!!

        $site = $chr . "\t" . $start_nt . "\t" . $stop_nt;
        $seen_elements{$name} = $site;
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't process input: $input\n";
    }
}
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

open my $BLASTN, '<', $blastn or die "Can't open BlastN file $blastn: $!";
while (my $input = <$BLASTN>) { 
    my $name = q{};

# Sample input:
# Query= ce6_ct_phastv3filtered_chrI.22 range=chrI:19498-19566 5'pad=0

    if ( $input =~ /\A Query = \s (\S+) /xms ) { 
        $name = $1;
        if (! exists $seen_elements{$name} ) { 
           die "Failed to get FASTA-header coordinates for element $name in BlastN file $blastn!\n";
        }
        # Record for the BlastN hits what they don't make easy to just read -- the sites:
        $blastn_sites{$seen_elements{$name}} = 1;
    }
    elsif ( $input =~ /\A Query = /xms ) { 
        die "Can't process input from BlastN file $blastn: $input\n";
    }
}
close $BLASTN or die "Can't close filehandle to BlastN file $blastn: $!";

open my $BED, '<', $bed or die "Can't open BED file $bed: $!";
while (my $input = <$BED>) { 
    chomp $input;
    my $site = q{};

# Sample input:
# chrIV	4930943	4931058	chrIV.8293	0	+

    if ( $input =~ /(chr \w+ \t \d+ \t \d+)/xms ) { 
        $site = $1;
    }
    if ( ( exists $blastn_sites{$site} ) and $badprint ) { 
        print "$input\n";
    }
    if ( $site and (! exists $blastn_sites{$site} ) and (! $badprint) ) {
        print "$input\n";
    }
}
close $BED or die "Can't close filehandle to BED file $bed: $!";

