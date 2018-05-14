#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();
my @mod_names   = ();

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
        # Note that '.*' flanks allow the OrderedLocusNames to be anywhere in the line,
        if ( $input =~ /\A GN \s+ .* OrderedLocusNames= (\S+) ;/xms ) {
            $mod_name = $1;
            # For TAIR, it will make life significantly easier if I switch the MOD IDs to full uppercase.
            $mod_name =~ tr/[a-z]/[A-Z]/;
            push @mod_names, $mod_name;
        }

        if ( $input =~ /\A GN \s+ Name= (\S+) ;/xms ) { 
            $gene_name = $1;
            push @gene_names, $gene_name;
        }
    }
}

# Fail-safe, but shouldn't be needed.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # We will only print data for UniProt proteins which have either an ENSEMBL or a human-readable gene name.
        if (@gene_names or @mod_names) {
            my @print_lines = ();

            # The description, at any rate, should be invariant per UniProt entry.
            $description =~ s/\A["]+//;
            $description =~ s/["]+\z//;
            if ( $description =~ /\S/xms ) {
                $description = "\"$description\"";
            }

            # If we do have gene names:
            if (@gene_names) { 
                foreach my $gene_name1 (@gene_names) { 
                    my $long_name = q{};
                    if (! @mod_names) { 
                        $long_name = $gene_name1;
                        my $print_text = "$uniprot\t$long_name\t$description";
                        push @print_lines, $print_text;
                    }
                    if (@mod_names) {
                        foreach my $mod_name1 (@mod_names) { 
                            $long_name = $mod_name1;
                            if ( $gene_name1 =~ /\S/xms ) { 
                                $long_name = $long_name . q{|} . $gene_name1;
                                my $print_text = "$uniprot\t$long_name\t$description";
                                push @print_lines, $print_text;
                            }
                        }
                    }
                }
            }

            # If we do not have gene names (but do, necessarily, have MOD ID names):
            if (! @gene_names) { 
                foreach my $mod_name1 (@mod_names) {
                    my $long_name = $mod_name1;
                    my $print_text = "$uniprot\t$long_name\t$description";
                    push @print_lines, $print_text;
                }
            }
            foreach my $print_text (@print_lines) {
                print "$print_text\n";
            }
        }
    }

    # With lines printed, return everything to zero for the next entry.
    @uniprots    = ();   
    @gene_names  = ();
    @mod_names   = ();
    
    $description = q{};   
    $gene_name   = q{};
    $mod_name    = q{};

    return;
}

