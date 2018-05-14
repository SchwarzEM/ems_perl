#!/usr/bin/env perl

# make_Csp7_Csp9_cDNA_lists_23aug2011.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/23/2011.
# Purpose: given two cd-hit-est-2d files of sp. 7 and sp. 9, find *all* sp. 7 and sp. 9 sequences showing cross-hits, and make them into sp. 7 and sp. 9 lists.

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $seq      = q{};
my @sp7_hits = ();
my @sp9_hits = ();

my @all_sp7_hits = ();
my @all_sp9_hits = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > Cluster \s+ \d+ \s* \z/xms ) { 
        if (@sp7_hits and @sp9_hits) {
            push @all_sp7_hits, @sp7_hits;
            push @all_sp9_hits, @sp9_hits;
        }
        @sp7_hits = ();
        @sp9_hits = ();
    }
    elsif ( $input =~ / > (Csp7.\S+_Length_\d+) /xms ) { 
        $seq = $1;
        push @sp7_hits, $seq;
    }
    elsif ( $input =~ / > (Csp9.\S+_Length_\d+) /xms ) {
        $seq = $1;
        push @sp9_hits, $seq;
    }
}

if (@sp7_hits and @sp9_hits) {
    push @all_sp7_hits, @sp7_hits;
    push @all_sp9_hits, @sp9_hits;
}

@all_sp7_hits = sort @all_sp7_hits;
@all_sp7_hits = uniq @all_sp7_hits;

@all_sp9_hits = sort @all_sp9_hits;
@all_sp9_hits = uniq @all_sp9_hits;

my $date = join('.', &get_local_date());

my $sp7_outfile = 'Csp7_hits_' . $date . '.txt';
my $sp9_outfile = 'Csp9_hits_' . $date . '.txt';

open my $SP7_LIST, '>', $sp7_outfile or die "Can't open sp. 7 list $sp7_outfile: $!";
foreach my $sp7_hit (@all_sp7_hits) { 
    print $SP7_LIST "$sp7_hit\n";
}
close $SP7_LIST or die "Can't close filehandle to sp. 7 list $sp7_outfile: $!";

open my $SP9_LIST, '>', $sp9_outfile or die "Can't open sp. 9 list $sp9_outfile: $!";
foreach my $sp9_hit (@all_sp9_hits) {
    print $SP9_LIST "$sp9_hit\n";
}
close $SP9_LIST or die "Can't close filehandle to sp. 9 list $sp9_outfile: $!";


sub get_local_date { 
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

