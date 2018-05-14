#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $subset_list = $ARGV[0];
my $full_table  = $ARGV[1];

my @seqs = ();
my %obs  = ();

if ( (! $subset_list) or (! $full_table) ) { 
    die "Format: extract_subset_repbase_table_24mar2014a.pl [subset list of desired sequences] [full data table of repeats] > [subset data table] ;\n";
}

open my $SUBSET, '<', $subset_list;
while (my $input = <$SUBSET>) { 
    chomp $input;
    if ( $input=~ /\A \S+ \z/xms ) { 
        push @seqs, $input;
        $obs{$input} = 1;
    }
    else {
        die "From subset_list file $subset_list, can't parse: $input\n";
    }
}
close $SUBSET;

open my $FULL, '<', $full_table;
while (my $input = <$FULL>) { 
    chomp $input;
    if ( $input=~ /\A (\S+) \t [^\t]+ \t [^\t]+ \t [^\t]+ \S \z/xms ) { 
        my $seqname = $1;
        if ( exists $obs{$seqname} ) { 
            print "$input\n";
        }
    }
    else { 
        die "From full table $full_table, can't parse: $input\n";
    }
}
close $FULL;

