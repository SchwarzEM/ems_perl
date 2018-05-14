#!/usr/bin/env perl

# bulyk_2_wb_pwm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/18/2009.
# Purpose: convert one or more Bulyk *.pwm files into a single rough-draft .ace file.

use strict;
use warnings;

my $i           = 78;    # HARD-CODED and totally contingent on matrices 1-77 already having been done.
my $wb_i        = q{};
my $protein     = q{};
my $residue     = q{};
my $raw_weights = q{};
my $rnd_weights = q{};
my @infiles = sort @ARGV;

foreach my $infile (@infiles) { 
    open my $INFILE, '<', $infile 
        or die "Can't open input file $infile!\n";
    print "\n";
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input =~ /\A Protein: \s+ (\S+) \s+ Seed \s+ k\-mer: \s+ [ACGT\.]+ \s+ Enrichment \s+ Score: \s+ \d \. \d+ \s* \z/xms ) { 
            $protein = $1;
            $wb_i = sprintf "%08u", $i;
            print "Position_Matrix : \"WBPmat$wb_i\"\n";
            if ( $protein =~ /\A ([\w\-]+)_([\w\-]+) \z /xms ) { 
                $protein = $1 . '/' . $2 . ' heterodimer';
            }
            print "Description       \"$protein binding site\.\" Paper_evidence \"WBPaper00034761\" // pmid19632181\n"; 

            # Note: Bulyk's 'pwm' files are actually PFMs!
            print "Type              Frequency\n";
        }
        elsif ( $input =~ / ([ACGT]) : ( (?: \s+ \d\.\d+ )+ ) \s* /xms ) { 
            $residue = $1;
            $raw_weights = $2;
            $rnd_weights =  $raw_weights;
            $rnd_weights =~ s/(\d\.\d{4})\d*/$1/g;
            $rnd_weights =~ s/\A\s+//;
            $rnd_weights =~ s/\s+\z//;
            $rnd_weights =~ s/(\s+)/  /g;
            print "Site_values       $residue  $rnd_weights\n";
        }
        elsif ( $input =~ /\S/xms ) { 
            die "Can't parse input line: $input\n";
        }
    }
    print "Remark            \"Original PWM $infile downloaded by curator from thebrain.bwh.harvard.edu/uniprobe/downloads/Cell09.\"\n";
    print "\n";
    $i++;
    close $INFILE or die "Can't close filehandle to input file $infile!\n";
}

