#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = '##gff-version 3';

while ( my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A[#]/xms ) {
        $input =~ s/\s+\z//;
        if ( $input =~ /\A ((?:\S+ \s+){28}) (.+) \z/xms ) {
            my $front = $1;
            my $desc  = $2;
            $front =~ s/\s+\z//;
            $desc  =~ s/\s+\z//;
            my @vals = split /\s+/, $front;
            my $rfam_motif  = $vals[1];
            my $rfam_acc    = $vals[2];
            my $seqname     = $vals[3];
            my $start_nt    = $vals[9];
            my $end_nt      = $vals[10];
            my $strand      = $vals[11];
            my $e_value     = $vals[17];

            print "$header\n" if $header;
            $header = q{};

            # Column 1: "seqid"
            print "$seqname";
            print "\t";

            # Column 2: "source"
            print "INFERNAL/Rfam";
            print "\t";

            # Column 3: "type"
            print "RNA_motif";
            print "\t";

            # Column 4: "start"
            print "$start_nt";
            print "\t";

            # Column 5: ""end"
            print "$end_nt";
            print "\t";

            # Column 6: "score"
            print "$e_value";
            print "\t";

            # Column 7: "strand"
            print "$strand";
            print "\t";

            # Column 8: "phase"
            print ".";
            print "\t";

            # Column 9: "attributes"
            print "ID=$rfam_acc|$rfam_motif; Name=\"$desc\"";
            print "\n";
        }
        else {
            die "Cannot parse: $input\n";
        }
    }
}
