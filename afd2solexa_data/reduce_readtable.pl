#!/usr/bin/env perl

# reduce_readtable.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/3/2011.
# Purpose: given one or more WS220-mapped read tables, which may have two entries per gene, pick entry with most total reads.

use strict;
use warnings;

my $data_ref;
my $gene         = q{};
my $data         = q{};
my $allreads     = q{};
my $length       = q{};
my @header_lines = ();

while (my $input = <>) { 
    chomp $input; 
    if ( $input !~ / \A WBGene\d+\S* \t \d+ \t \d+ \t \d+ \t \d+ \t \d+ \t \d+\.\d+ \z /xms ) { 
        push @header_lines, $input;
    }
    if ( $input =~ / \A (WBGene\d+\S*) \t ((\d+) \t \d+ \t \d+ \t \d+ \t \d+ \t (\d+\.\d+) ) \z /xms ) {
        $gene     = $1;
        $data     = $2;
        $allreads = $3;
        $length   = $4;
        if (    (! $data_ref->{$gene}->{'allreads'}                     ) 
             or ( $allreads > $data_ref->{$gene}->{'allreads'}          ) 
             or (     ( $allreads == $data_ref->{$gene}->{'allreads'} ) 
                  and ( $length  > $data_ref->{$gene}->{'length'}     ) ) ) { 
            $data_ref->{$gene}->{'allreads'} = $allreads;
            $data_ref->{$gene}->{'length'}   = $length;
            $data_ref->{$gene}->{'data'}     = $data;
        }
    }
}

# Warn if there are multiple lines getting treated as header lines.
my $header_line_count = @header_lines;
if ( $header_line_count > 1 ) { 
    warn "Multiple header lines: @header_lines\n";
}

my $final_header_line = $header_lines[0];

# The following is a kludge to deal with headers from an earlier step of the pipeline that I can't rerun cleanly on 6/3/2011.
# I commented it out on 7/1/2012, because it's awkward -- I don't want or need the "Gene" header for DESeq, and I can use
#     in-line Perl edits if I do need it for other things anyway.
#
# if ( $final_header_line =~ /\A \t/xms) { 
#     $final_header_line = 'Gene' . $final_header_line;
# }

print "$final_header_line\n";
foreach my $wbgene ( sort keys %{ $data_ref } ) { 
    print "$wbgene\t$data_ref->{$wbgene}->{'data'}\n";
}

