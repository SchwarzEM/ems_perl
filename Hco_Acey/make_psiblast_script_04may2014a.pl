#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq);

my $query_file = $ARGV[0];
my $e_value    = $ARGV[1];
my $iterations = $ARGV[2];

if ( (! $query_file) or (! $e_value) or (! $iterations) ) { 
    die "Format: make_psiblast_script_04may2014a.pl [file with query list] [E-value threshold] [iteration number]\n";
}

if ( (! $e_value) or ( $e_value <= 0 ) or (! looks_like_number($e_value) ) ) { 
    die "E-value needs to be positive real number, not this: $e_value\n";
}

if ( ( $iterations != int($iterations) ) or ( $iterations <= 0 ) ) {
    die "Iteration number needs to be positive integer, not this: $iterations\n";
}

my @queries = qw();

open my $QUERY_FILE, '<', $query_file;
while ( my $query = <$QUERY_FILE> ) {
    chomp $query;
    $query =~ s/\A\s+//;
    $query =~ s/\s+\z//;
    if ( $query =~ /\s/xms ) { 
        die "In query list file $query_file, can't parse query: $query\n";
    }
    push @queries, $query;
}

@queries = uniq(@queries);

my $header = '#!/bin/bash' . "\n\n";
my $footer = "    echo done_psiblast_script_04may2014a > script_04may2014a.txt ;\n\n";

foreach my $query (@queries) { 
    print $header if $header;
    $header = q{};

    my $stem = $query;
    $stem = basename($stem);
    $stem    =~ s/\.fa\z//;

    my $output1 = "$stem.psiblast_" . "$iterations" . 'x_' . "$e_value.metazoan_prots.tsv.txt";
    my $output2 = "$stem.psiblast_" . "$iterations" . 'x_' . "$e_value.metazoan_prots.final.tsv.txt";
    my $output3 = "$stem.psiblast_" . "$iterations" . 'x_' . "$e_value.metazoan_prots.final.hits.txt";

    print "    psiblast -db ../dbs/metazoan_prots_14apr2014",
          " -query $query -num_threads 8 -evalue $e_value -inclusion_ethresh $e_value -num_iterations $iterations -outfmt 6",
          " -seg yes -out $output1 ;\n",
          ;
    print "    cat $output1 | get_final_psiblast_tsv.pl > $output2 ;\n";
    print "    cut -f 2 $output2 | sort | uniq > $output3 ;\n";
    print "\n";
}
print $footer if (! $header);

