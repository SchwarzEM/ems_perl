#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my @uniprots    = ();
my @gene_names  = ();

my %hum2cds    = ();
my %cds2hum    = ();

my $description = q{};
my $gene_name   = q{};

my $cds_name    = q{};
my $human_name  = q{};

# Sample input -- but note that there can be two or more different GN, etc. per protein.
# 
# ID   2A5D_YEAST              Reviewed;         757 AA.
# AC   P38903; D6W281;
# [...]
# DE   RecName: Full=Serine/threonine-protein phosphatase 2A 56 kDa regulatory subunit delta isoform;
# [...]
# [various different instances:]
# GN   Name=CDC55; OrderedLocusNames=YGL190C; ORFNames=G1345;
# GN   OrderedLocusNames=YJR149W; ORFNames=J2213;
# GN   Name=BNA1; Synonyms=HAD1; OrderedLocusNames=YJR025C; ORFNames=J1550;

# One tough input:
# ID   AFI1_YEAST              Reviewed;         893 AA.
# AC   Q99222; D6W2I7; Q92274;
# GN   Name=AFI1; OrderedLocusNames=YOR129C; ORFNames=YOR3296C;
# DR   SGD; S000005655; YOR129C.

# Yet another:
# GN   OrderedLocusNames=YJR149W; ORFNames=J2213;

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
    # Note that '.*' flanks allow the OrderedLocusNames to be anywhere in the line.
    elsif ( @uniprots and ( $input =~ /\A GN \s+ /xms ) ) { 
        # Sometimes happens on the same line as the CDS, other times one line ahead...
        if ( $input =~ /\A GN \s+ Name= (\S+) ;/xms ) {
            $human_name = $1;
            # Deal with inconsistencies of lower-case here vs. upper-case down in the SGD line:
            $human_name =~ tr/[a-z]/[A-Z]/;
        }

        # This is the real CDS data.  Assume that it will never happen *before* we see the human-readable gene.
        if ( $input =~ /\A GN \s+ .* OrderedLocusNames= (\S+) ;/xms ) { 
            $cds_name   = $1;

            # Unless we got the human name one line before, which sometimes stupidly happens...
            if (! $human_name) { 
                $human_name = $cds_name;
            }

            # Note that, later on, we will *mostly* need to interpret a human name next to an SGD number.
            $hum2cds{$human_name} = $cds_name;

            # But there are backwards cases where the SGD number comes with a CDS instead, so prepare to deal with them too.
            $cds2hum{$cds_name} = $human_name;
        }
    }
    elsif ( @uniprots and ( $input =~ /\A DR \s+ SGD; /xms ) ) { 
        if ( $input =~ /\A DR \s+ SGD; \s+ (S\d+) ; \s+ (\S+) \. \z /xms ) { 
            my $mod_name            = $1;
            my $putative_human_name = $2;
            my $putative_cds_name   = $putative_human_name;
            $gene_name     = q{};

            # There are, alas, many reasons this can happen.
            # I have resorted to keeping older values of $cds_gene and $human_gene around, as last resorts.
            if (! $hum2cds{$putative_human_name} ) { 
                # If I catch an error in ordering of CDS and human genes, try to switch things around.
                if ( $cds2hum{$putative_human_name} ) { 
                    my $real_human_name = $cds2hum{$putative_human_name};
                    $cds_name           = $putative_human_name;
                    $human_name         = $real_human_name;
                }
                # Otherwise, by default, fall back on the previously stored values for $cds_gene and $human_gene.
                $gene_name = $mod_name . q{|} . $cds_name . q{|} . $human_name;
            }

            else { 
                if ( $putative_human_name eq $hum2cds{$putative_human_name} ) {
                    $gene_name = $mod_name . q{|} . $putative_cds_name;
                }
                if ( $human_name ne $hum2cds{$putative_human_name} ) {
                    $putative_cds_name = $hum2cds{$putative_human_name};
                    $gene_name         = $mod_name . q{|} . $putative_cds_name . q{|} . $putative_human_name;
                }
            }
            push @gene_names, $gene_name;
        }  
        else {
            die "Can't parse input: $input\n"; 
        }
    }
}

# Fail-safe, but shouldn't be needed.
map_print_and_clear_stored_data();

# This has been recoded to explicitly handle any possible mixture of 0 to 2+ gene names or ENSEMBL names per protein.
sub map_print_and_clear_stored_data {
    foreach my $uniprot (@uniprots) { 
        # Only @gene_names which are already as long as possible, so life is simple.
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

    %hum2cds     = ();
    %cds2hum     = ();

    $cds_name    = q{};
    $human_name  = q{};

    return;
}

