#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();

my $description = q{};
my $gene_name   = q{};
my $mod_name    = q{};

# Sample input -- but note that there can be two or more different GN, etc. per protein.
# Also note that UniProt does not always agree with TAIR about how to name a gene (e.g., LGO).
# 
# ID   SMR1_ARATH              Reviewed;         128 AA.
# AC   Q9LPP4; Q8GXB8; Q8LDA2;
# [...]
# DE   RecName: Full=Cyclin-dependent protein kinase inhibitor SMR1;
# [...]
# GN   Name=SMR1; Synonyms=LGO; OrderedLocusNames=At3g10525;
# GN   ORFNames=F13M14.31, F18K10.10;
# [......]
# ID   SMR3_ARATH              Reviewed;         115 AA.
# AC   Q9LZ60; Q8LD37;
# [...]
# DE   RecName: Full=Cyclin-dependent protein kinase inhibitor SMR3;
# [...]
# GN   Name=SMR3; OrderedLocusNames=At5g02420; ORFNames=T22P11.10;

# In contrast, my mapping from TAIR data gave me:
# Q9LPP4	AT3G10525|LGO

# The thing to do, probably, is let UniProt names stand while I am doing the extraction,
#     but then 'polish' the results to fit TAIR standards later.  Sigh, more names.

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\s+\z//;

    # '//' == end of entry.
    if ( $input =~ /\A \/\/ /xms ) { 
        if (@uniprots) {
            map_print_and_clear_stored_data();
        }
    }

    # Note that there can be more than one AC line per entry; important not to overwrite first line's accessions with the 2+ lines.
    if ( $input =~ /\A AC \s+ (\S.*\S) \z/xms ) { 
        my $uniprot_text = $1;
        $uniprot_text =~ s/;//g;
        my @new_uniprots = split /\s+/, $uniprot_text;
        push @uniprots, @new_uniprots;
    }
    elsif ( @uniprots and ( $input =~ /\A DE \s+ RecName: \s+ Full= (.+) ; \z/xms ) ) {
        $description = $1;
    }
    elsif ( @uniprots and ( $input =~ /\A GN \s+ /xms ) ) { 
        # Note that '.*' flanks allow the OrderedLocusNames to be anywhere in the line.
        # For TAIR, it will make life significantly easier if I switch all gene names or gene name components to full uppercase.
        if ( $input =~ /\A GN \s+ Name= (\S+) ; .* OrderedLocusNames= (\S+) ;/xms ) {
            $gene_name = $1;
            $mod_name  = $2;
            # For TAIR, it will make life significantly easier if I switch the MOD IDs to full uppercase.
            $gene_name = $mod_name . q{|} . $gene_name;
            $gene_name =~ tr/[a-z]/[A-Z]/;
            push @gene_names, $gene_name;
        }
        elsif ( $input =~ /\A GN \s+ .* OrderedLocusNames= (\S+) ;/xms ) {
            $mod_name = $1;
            # For TAIR, it will make life significantly easier if I switch the MOD IDs to full uppercase.
            $mod_name =~ tr/[a-z]/[A-Z]/;
            push @gene_names, $mod_name;
        }
        elsif ( $input =~ /\A GN \s+ Name= (\S+) ;/xms ) {
            $gene_name = $1;
            $gene_name =~ tr/[a-z]/[A-Z]/;
            push @gene_names, $gene_name;
        }
    }
}

# Fail-safe, but shouldn't be needed.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # Keeping life simple, we have made all gene names long names, and do not separately track mod names.
        if (@gene_names) {
            # The description, at any rate, should be invariant per UniProt entry.
            $description =~ s/\A["]+//;
            $description =~ s/["]+\z//;
            if ( $description =~ /\S/xms ) {
                $description = "\"$description\"";
            }

            foreach my $gene_name1 (@gene_names) { 
                print "$uniprot\t$gene_name1\t$description\n";
            }
        }
    }

    # With lines printed, return everything to zero for the next entry.
    @uniprots    = ();   
    @gene_names  = ();
    
    $description = q{};   
    $gene_name   = q{};
    $mod_name    = q{};

    return;
}

