#!/usr/bin/env perl

# sort_dge_groups_06feb2016.pl -- Erich Schwarz <ems394@cornell.edu>, 2/6/2016.
# Purpose: given gene annot file with well-headered fold-changes and padj, automatically parse and export sig. DEG groups.

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $data_ref;

while (my $input = <>) {
    chomp $input;
    # from the header line, identify which columns contain which data types.
    if ( $input =~ /\A (\S+) \t (.+\S) \z/xms ) { 
        my $gene_id    = $1;
        my $field_text = $2;
        my $i          = 0;
        my @fields = split /\t/, $field_text;
        foreach my $field (@fields) {
            $i++;
            if ( ( $gene_id eq 'Gene' ) and ( $field =~ /\A (log2FoldChange|padj) \[ (\S+) \] \z/xms ) ) {
                my $data_type = $1;
                my $condition = $2;
                $data_ref->{'field'}->{$i}->{'data_type'} = $data_type;
                $data_ref->{'field'}->{$i}->{'condition'} = $condition;
            }
            elsif ( $gene_id ne 'Gene' ) { 
                if (     ( exists $data_ref->{'field'}->{$i}->{'data_type'} ) 
                     and ( exists $data_ref->{'field'}->{$i}->{'condition'} ) ) { 

                    my $data_type = $data_ref->{'field'}->{$i}->{'data_type'};                
                    my $condition = $data_ref->{'field'}->{$i}->{'condition'};

                    # For genes w/o sig. changes in given comparison, may well be no valid data in the field; 
                    #     no point in doing the mapping unless there is!
                    if (looks_like_number($field) ) {
                        # Order the data in a way that makes it easy to mine efficiently later:
                        $data_ref->{'condition'}->{$condition}->{'gene'}->{$gene_id}->{'data_type'}->{$data_type} = $field;
                    }
                }
            }
        }
    } 
    else { 
        die "Cannot parse input line: $input\n";
    }
}

my @interesting_conds = sort keys %{ $data_ref->{'condition'} };

foreach my $int_cond (@interesting_conds) {
    my $upfile    = "$int_cond.upreg.tsv.txt";
    my $downfile  = "$int_cond.downreg.tsv.txt";
    $upfile       = safename($upfile);
    $downfile     = safename($downfile);

    my $up_header = "Gene\tlog2FoldChange\tpadj";
    my $dn_header = $up_header;

    open my $UP,   '>', $upfile;
    open my $DOWN, '>', $downfile;

    my @genes = sort keys %{ $data_ref->{'condition'}->{$int_cond}->{'gene'} };

    foreach my $gene (@genes) {
        my $log2FoldChange = $data_ref->{'condition'}->{$int_cond}->{'gene'}->{$gene}->{'data_type'}->{'log2FoldChange'};
        my $padj = $data_ref->{'condition'}->{$int_cond}->{'gene'}->{$gene}->{'data_type'}->{'padj'};
        if ( $log2FoldChange > 0 ) {
            print $UP "$up_header\n" if $up_header;
            $up_header = q{};
            print $UP "$gene\t$log2FoldChange\t$padj\n";
        }
        if ( $log2FoldChange < 0 ) {
            print $DOWN "$dn_header\n" if $dn_header;
            $dn_header = q{};
            print $DOWN "$gene\t$log2FoldChange\t$padj\n";
        }
    }
    close $UP;
    close $DOWN;
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

