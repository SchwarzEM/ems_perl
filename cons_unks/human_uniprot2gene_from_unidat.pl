#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();
my @ens_names   = ();

my $description = q{};
my $gene_name   = q{};
my $ens_name    = q{};

# Sample input -- but note that there can be two or more different GN and ENSG per protein!!
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

    # '//' == end of record.
    if ( $input =~ /\A \/ \/ /xms ) { 
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
    elsif ( @uniprots and ( $input =~ /\A GN \s+ Name= (\S+) ;/xms ) ) { 
        $gene_name = $1;
        push @gene_names, $gene_name;
    }

    # Note that '.*' flanks allow the (ENSG\d+) to be anywhere in the line, which it can easily be.
    elsif ( @uniprots and ( $input =~ /\A DR \s+ Ensembl; \s+ .* (ENSG\d+) .* \z/xms ) ) { 
        $ens_name = $1;
        push @ens_names, $ens_name;
    }
}

# Clear out last stored data, after the text is done running.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # We will only print data for UniProt proteins which have either an ENSEMBL or a human-readable gene name.
        if (@gene_names or @ens_names) {
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
                    if (! @ens_names) { 
                        $long_name = $gene_name1;
                        my $print_text = "$uniprot\t$long_name\t$description";
                        push @print_lines, $print_text;
                    }
                    if (@ens_names) {
                        foreach my $ens_name1 (@ens_names) { 
                            $long_name = $ens_name1;
                            if ( $gene_name1 =~ /\S/xms ) { 
                                $long_name = $long_name . q{|} . $gene_name1;
                                my $print_text = "$uniprot\t$long_name\t$description";
                                push @print_lines, $print_text;
                            }
                        }
                    }
                }
            }

            # If we do not have gene names (but do, necessarily, have ENSEMBL names):
            if (! @gene_names) { 
                foreach my $ens_name1 (@ens_names) {
                    my $long_name = $ens_name1;
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
    @ens_names   = ();
    
    $description = q{};   
    $gene_name   = q{};
    $ens_name    = q{};

    return;
}

