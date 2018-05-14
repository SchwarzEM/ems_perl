#!/usr/bin/env perl

use strict;
use warnings;

my $motif = $ARGV[0];
my $table = $ARGV[1];
my $prots = $ARGV[2];

my %ok_prots = ();

my $header = '#' . " Sequence\tMatch_length\tMotif_length [$motif]\tStart_aa\tStop_aa";

if (! $motif) {
    die "Format: extract_pfam_matchzones_18nov2013.pl [motif] [domain table] [list of OK proteins]\n";
}

if ( $motif !~ /\A \S+ \z/xms ) { 
    die "Need single word for motif, not $motif\n";
}

open my $PROTS, '<', $prots or die "Can't open protein list $prots: $!";
while (my $input = <$PROTS>) {
    chomp $input;
    if ( $input !~ /\A \S+ \z/xms ) {
        die "From protein list $prots, cannot parse: $input\n";
    }
    $ok_prots{$input} = 1;
}
close $PROTS or die "Can't close filehandle to protein list $prots: $!";

open my $TABLE, '<', $table or die "Can't open PFAM domain table $table: $!";
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ \S+ \s+ \d+ \s+ \S+ \s+ (\S+) \s+ (\d+) .+ \s+ \d+ \s+ \d+ \s+ (\d+) \s+ (\d+) \s+ \d+ \s+ \d+ \s+ \d\.\d{2} \s+ \S+ \s* \z/xms ) { 
        my $prot      = $1;
        my $obs_motif = $2;
        my $motif_len = $3;
        my $start_aa  = $4;
        my $end_aa    = $5;
        if ( $ok_prots{$prot} and ( $motif eq $obs_motif ) ) {
            my $size_of_match = ($end_aa - $start_aa);
            # Only print header once at the top, *if* there is something else to print.
            if ($header) {
                print "$header\n";
                $header = q{};
            }
            print "$prot\t$size_of_match\t$motif_len\t$start_aa\t$end_aa\n";
        }
    }
}
close $TABLE or die "Can't close filehandle to PFAM domain table $table: $!";

