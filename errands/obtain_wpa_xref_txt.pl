#!/usr/bin/perl

# obtain_wpa_xref_txt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/13/2007.
# Purpose: extract handy reference table in ASCII from somewhat klunky wpa_xref.cgi site.

use strict;
use warnings;
use LWP::Simple;

my $wbref_url = 'http://tazendra.caltech.edu/~postgres/cgi-bin/wpa_xref.cgi';
my $wbref_list = get($wbref_url);

unless (defined $wbref_list) {
    die "Failed to download $wbref_url\n";
}

my @wbrefs = split /\n/, $wbref_list;

my %wb_papers    = ();
my %cgc_aliases  = ();
my %pmid_aliases = ();

foreach my $wbref (@wbrefs) { 
    if ($wbref =~ /(WBPaper\d+).+(cgc\d+)/) {
        $wb_papers{$1} = 1;
        $cgc_aliases{$1} = $2;
    }
    elsif ($wbref =~ /(WBPaper\d+).+(pmid\d+)/) {
        $wb_papers{$1} = 1;
        $pmid_aliases{$1} = $2;
    }
}

foreach my $wb_paper (sort keys %wb_papers) { 
    print "[$wb_paper";
    if ($cgc_aliases{$wb_paper}) {
        print "; $cgc_aliases{$wb_paper}";
    }
    if ($pmid_aliases{$wb_paper}) {
        print "; $pmid_aliases{$wb_paper}";
    }
    print "]\n";
}

