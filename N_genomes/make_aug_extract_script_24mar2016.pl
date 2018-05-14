#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '#!/bin/bash' . "\n\n";

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $seq  = $1;
        my $gff  = $2;
        my $stem = q{};

        if ( (! -r $seq) or (! -r $gff ) ) {
            die "Cannot read at least one of these files: $input\n";
        }

        if ( $gff =~ /\A (\S+) \.gff \z/xms ) { 
            $stem = $1;
        }
        else {
            die "Cannot extract stem from GFF name.\n";
        }

        print $header if $header;
        $header = q{};

        print '    /mnt/home/emsch/src/augustus-3.2.1/scripts/getAnnoFasta.pl --seqfile=';
        print "$seq";
        print " $gff ;\n";
        print "\n";

        print "    uniform_fasta.pl -a -k -i $stem.aa > $stem.pep.fa ;\n";
        print "    uniform_fasta.pl -a -k -i $stem.codingseq > $stem.cds_dna.fa ;\n";
        print "    uniform_fasta.pl -a -k -i $stem.cdsexons > $stem.exons_dna.fa ;\n";
        print "\n";
        
        print "    rm $stem.aa $stem.codingseq $stem.cdsexons ;\n";
        print "\n";

        print "    get_largest_isoforms.pl -t aug -i $stem.pep.fa > $stem.max_iso.pep.fa ;\n";
        print "    get_largest_isoforms.pl -t aug -i $stem.cds_dna.fa > $stem.max_iso.cds_dna.fa ;\n";
        print "\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
}

