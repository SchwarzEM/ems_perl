#!/usr/bin/env perl

# fastq_adapt_trim.pl -- Ross Hall <rossh@unimelb.edu.au>, 7/1/2010.
# Purpose: trims quality of FASTQ file reads, apparently.
#
# Detailed purpose:
#
# Script reads sequences from a FASTQ file and
# searches for location of motif.
# If motif found > adaptor_threshold_location print all bases up to motif
# If motif found < adaptor_threshold_location ignore sequence
#
# Trims quality to same as sequence
# If motif not found - trim to max_length bases wide
#
# Usage: fastq_trim.pl motif adaptor_threshold_location max_length input_FASTQ_file > output_faq_file

use strict;
use warnings;

use lib '/usr/share/perl5/';
use Bio::Perl;
use Bio::SeqIO;
use Bio::DB::GenBank;

if ( @ARGV != 4 ) {
    print STDERR
"fasta_trim.pl motif adaptor_threshold_location max_length input_FASTA_file\n";
    exit(1);
}

my $motif     = shift;
my $thresh    = shift;
my $maxlength = shift;
my $infile    = shift;

open( INFILE, $infile ) || &ErrorMessage( "Cannot open file " . $infile );

my $line;
my $motif_pos = 0;
my $sub;
my $inseq  = 0;
my $inqual = 0;
my $seqheader;
my $qualheader;
my $printseq = 0;

while (<INFILE>) {
    $line = $_;
    chomp($line);
    if ( $line =~ /^@/ && $line =~ /[0-9]/ ) {

        # sequence header
        $seqheader = $line;
        $inseq     = 1;
        next;
    }
    if ( $inseq == 1 ) {

        # sequence line
        $motif_pos = length($line);
        $printseq  = 1;
        if ( $line =~ /$motif/ ) {

            # motif is in the sequence
            # get its location
            $motif_pos = index( $line, $motif );
            if ( $motif_pos < $thresh ) {
                $printseq = 0;
            }
        }

        if ($printseq) {
            if ( $motif_pos > $maxlength ) { $motif_pos = $maxlength; }
            $sub = substr( $line, 0, $motif_pos );
            print "$seqheader\n";
            print "$sub\n";
        }
        $inseq = 0;
        next;
    }

    if ( $line =~ /^\+/ ) {

        # quality header
        $qualheader = $line;
        $inqual     = 1;
        next;
    }

    if ( $inqual == 1 ) {

        # quality line - same length as sequence
        if ($printseq) {
            if ( $motif_pos > $maxlength ) { $motif_pos = $maxlength; }
            $sub = substr( $line, 0, $motif_pos );
            print "$qualheader\n";
            print "$sub\n";
        }
        $inqual = 0;
    }
}

close(INFILE);

sub ErrorMessage {
    my $msg = shift;
    print "Fatal error: $msg\n";
    exit(1);
}

