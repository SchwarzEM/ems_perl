#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my $description = q{};
my $gene_name   = q{};
my $ens_name    = q{};

# Sample input:
# 
# ID   1433B_HUMAN             Reviewed;         246 AA.
# AC   P31946; A8K9K2; E1P616;
# [...]
# DE   RecName: Full=14-3-3 protein beta/alpha;
# [...]
# GN   Name=YWHAB;
# [...]
# DR   Ensembl; ENST00000353703; ENSP00000300161; ENSG00000166913.

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s+\z//;
    if ( $input =~ /\A ID \s+ \S+ /xms ) { 
        if (@uniprots) {
            map_print_and_clear_stored_data();
        }
    }
    if ( $input =~ /\A AC \s+ (\S.*\S) \z/xms ) { 
        my $uniprot_text = $1;
        $uniprot_text =~ s/;//g;
        @uniprots = split /\s+/, $uniprot_text;
    }
    elsif ( @uniprots and ( $input =~ /\A DE \s+ RecName: \s+ Full= (.+) ; \z/xms ) ) {
        $description = $1;
    }
    elsif ( @uniprots and ( $input =~ /\A GN \s+ Name= (\S+) ;/xms ) ) { 
        $gene_name = $1;
    }
    elsif ( @uniprots and ( $input =~ /\A DR \s+ Ensembl; \s+ .* (ENSG\d+) .* \z/xms ) ) { 
        $ens_name = $1;
    }
}

# Clear out last stored data, after the text is done running.
map_print_and_clear_stored_data();

sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # We will only print data for UniProt proteins which have either an ENSEMBL or a human-readable gene name.
        if ($gene_name or $ens_name) {
            my $long_name = q{};
            if (! $ens_name) { 
                warn "Failed to map uniprot $uniprot to an ENSEMBL gene name\n";
                $long_name = $gene_name;
            }
            if ($ens_name) {
                $long_name = $ens_name;
                if ( $gene_name =~ /\S/xms ) { 
                    $long_name = $long_name . q{|} . $gene_name;
                }
            }
            if (! $description) {
                warn "Failed to map protein description for uniprot $uniprot (encoded by gene $long_name)\n";
            }
            print "$uniprot\t$long_name\t\"$description\"\n";
        }
    }
    @uniprots    = ();
    $description = q{};
    $gene_name   = q{};
    $ens_name    = q{};    
    return;
}

