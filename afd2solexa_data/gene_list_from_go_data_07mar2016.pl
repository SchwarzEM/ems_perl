#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %obs_genes = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \s \[(GO:\d+)\] \t (?: [^\t]* \t){4} ([^\t]+) /xms ) { 
        my $go_text   = $1;
        my $go_id     = $2;
        my $gene_text = $3;

        %obs_genes = ();

        $go_text =~ s/\A\s+//;
        $go_text =~ s/\s+\z//;
        $go_text =~ s/\s/_/g;
        $go_text =~ s/__/_/g;
        $go_text =~ s/([,]|[']|["]|[(]|[)]|[:])/./g;
        $go_text =~ s/[.]+/./g;

        $go_id   =~ s/[:]/_/g;

        my $go_gene_list = "$go_text.$go_id.gene_list.txt";
        $go_gene_list    = safename($go_gene_list);

        while ( $gene_text =~ /(AT (?: 1|2|3|4|5|C|M) G\d+)/xmsg ) {
            my $tair_gene = $1;
            $obs_genes{$tair_gene} = 1;
        }
        
        my @genes = sort keys %obs_genes;

        open my $LIST, '>', $go_gene_list;
        foreach my $gene (@genes) {
            print $LIST "$gene\n";
        }
        close $LIST;
    }
}


sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

