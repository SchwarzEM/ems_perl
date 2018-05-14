#!/usr/bin/env perl

# pfam2annot.pl -- Erich Schwarz <emsch@caltech.edu>, 11/9/2012.
# Purpose: given a gene2tx table and a PFAM 3.0 output, make a 2-column table of genes and PFAM domains.

use strict;
use warnings;
use Getopt::Long;

my $data_ref;

my $pfam    = q{};
my $help;

GetOptions ( 'pfam=s'    => \$pfam,
             'help'      => \$help, );

if ( $help or (! $pfam) ) { 
    die "Format: pfam_to_pfam_list.pl --pfam|-p [PFAM 3.0 TSV output] --help|-h [print this message]\n";
}

open my $PFAM, '<', $pfam or die "Can't open PFAM table: $pfam\n";
while (my $input = <$PFAM>) {
    chomp $input;
    # Acey_2012.08.05_0231.g3009.t2 -          1-cysPrx_C           PF10417.4 
    if ( $input !~ /\A [#] /xms ) { 
        if ( $input =~ /\A \S+ \s+ \S+ \s+ (\S+) \s+ (\S+) \s /xms ) { 
            my $pfam_name = $1;
            my $pfam_acc  = $2;
            my $pfam_desc = "$pfam_name [$pfam_acc]";
            $data_ref->{'pfam_desc'}->{$pfam_desc} = 1;
        }
        else { 
            die "From PFAM table $pfam, can't parse input: $input\n";
        } 
    }
}
close $PFAM or die "Can't close filehandle to PFAM table: $pfam\n";

my @pfam_descs = sort keys %{ $data_ref->{'pfam_desc'} };

foreach my $pfam_desc1 (@pfam_descs) {
    print "$pfam_desc1\n";
}

