#!/usr/bin/env perl

# fam2annots_14nov2018.pl -- Erich Schwarz <ems394@cornell.edu>, 11/14/2018.
# Purpose: given multiple gene2pmid and gene2[family] files, along with a key taxon and other required taxa, produce ranked tables of families with summed and per-gene annotation densities.

use strict;
use warnings;
use autodie;

use Getopt::Long;
use List::MoreUtils qw(uniq);

my @gene2pmids = ();
my @gene2fams  = ();
my $key_taxon  = q{};
my @req_taxa   = ();

my $header = "Gene"
             . "\tFamily"
             . "\tGene_summary"
             . "\tGene_list"
             . "\tGene_total"
             . "\tAnnot_sum"
             . "\tAnnot_sum/gene_total"
             ;

my $data_ref;

my $help;

GetOptions ( 'pmids=s{,}'    => \@gene2pmids,
             'fams=s{,}'     => \@gene2fams,
             'key=s',        => \$key_taxon,
             'req_taxa=s{,}' => \@req_taxa,
             'help'          => \$help,   );

if ( $help or (! @gene2pmids) or (! @gene2fams) or (! $key_taxon) ) { 
    die "Format: fam2annots_17nov2018.pl\n",
        "    --pmids|-p     [gene2pmid files]\n",
        "    --fams|-f      [gene2fam files]\n",
        "    --key|-k       [key taxon name]\n",
        "    --req_taxa|-r  [other required taxa]\n",
        "    --help|-h      [print this message]\n",
        ;
}

# Sample input line:
# elegans|WBGene00000490|che-11	10778742 [0.000101677681748856]; (...) 9851916 [3.81184722116338e-05]	1.71009965647908

# Store gene-to-pmid data.
foreach my $gene2pmid (@gene2pmids) {
    open my $PMID, '<', $gene2pmid;
    while (my $input = <$PMID>) {
        chomp $input;
        if ( $input =~ /\A (([a-z]+)\|\S+) \t ([^\t]+) \t (\S+) \z/xms ) {
            my $gene  = $1;
            my $taxon = $2;
            my $pmids = $3;
            my $annot = $4;

            if ( exists $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene} ) {
                die "Redundant mapping of taxon $taxon to gene $gene\n";
            }
            $data_ref->{'taxon'}->{$taxon}->{'gene'}->{$gene} = 1;

            if ( exists $data_ref->{'gene'}->{$gene}->{'pmids'} ) {
                die "For gene $gene, redundant PMIDs: \"$data_ref->{'gene'}->{$gene}->{'pmids'}\" versus \"$pmids\"\n";
            }
            $data_ref->{'gene'}->{$gene}->{'pmids'} = $pmids;

       	    if ( exists	$data_ref->{'gene'}->{$gene}->{'annot'} ) {
       	       	die "For gene $gene, redundant annot. weights: $annot versus $data_ref->{'gene'}->{$gene}->{'annot'}\n";
            }
            $data_ref->{'gene'}->{$gene}->{'annot'} = $annot;
        }
        else {
            die "In gene2pmid file $gene2pmid, cannot parse: $input\n";
        }
    }
    close $PMID;
}

# Sample input:
# elegans|WBGene00000001|aap-1	pfam|PF00017|SH2

# Store gene-to-family data.
foreach my $gene2fam (@gene2fams) {
    open my $FAM, '<', $gene2fam;
    while (my $input = <$FAM>) {
        chomp $input;
        if ( $input =~ /\A (([a-z]+)\|\S+) \t (\S+) \z/xms ) {
            my $gene   = $1;
            my $taxon  = $2;
            my $family = $3;

            # This will often have to happen multiple times, so tolerate that:
            $data_ref->{'family'}->{$family}->{'taxon'}->{$taxon} = 1;

            if ( exists $data_ref->{'gene'}->{$gene}->{'family'}->{$family} ) {
                die "Redundant mapping of gene $gene to family $family\n";
            }
            $data_ref->{'gene'}->{$gene}->{'family'}->{$family} = 1;

            if ( exists $data_ref->{'family'}->{$family}->{'gene'}->{$gene} ) {
                die "Redundant mapping of family $family to gene $gene\n";
            }
            $data_ref->{'family'}->{$family}->{'gene'}->{$gene} = 1;
        }
        else {
            die "In gene2fam file $gene2fam, cannot parse: $input\n";
        }
    }
    close $FAM;
}

# List genes that belong to $key_taxon and were mapped to at least one gene family.
my @key_genes = sort keys %{ $data_ref->{'taxon'}->{$key_taxon}->{'gene'} };

# Nonredundantly list gene families that contain a $key_taxon genes, and (optionally) have members from any other required taxa.
foreach my $key_gene (@key_genes) {
    my @init_fams = sort keys %{ $data_ref->{'gene'}->{$key_gene}->{'family'} };

    # Nonredundantly record all the families containing 1+ key genes.
    foreach my $init_family (@init_fams) {
        $data_ref->{'key_gene'}->{$key_gene}->{'key_family'}->{$init_family} = 1;
    }

    # *Optionally*, delete all the families that fail to contain at least one other @req_taxa.
    # If we leave @req_taxa empty, it doesn't happen.
    # Only delete a record if it still exists; don't try to delete 2+ times!
    foreach my $init_family (@init_fams) {
        if (@req_taxa) {
            foreach my $req_taxon (@req_taxa) {
                if (     ( exists $data_ref->{'key_gene'}->{$key_gene}->{'key_family'}->{$init_family} ) 
                     and (! exists $data_ref->{'family'}->{$init_family}->{'taxon'}->{$req_taxon}      ) ) {
                    delete $data_ref->{'key_gene'}->{$key_gene}->{'key_family'}->{$init_family};
                }
            }
        }
    }
}

# Redefine @key_genes so that they are required to retain at least one key_family.
@key_genes = grep { ( exists $data_ref->{'key_gene'}->{$_}->{'key_family'} ) } sort keys %{ $data_ref->{'key_gene'} };

# Compute the summed gene count and annot value for each family of interest.
foreach my $key_gene (@key_genes) {
    my @key_families = sort keys %{ $data_ref->{'key_gene'}->{$key_gene}->{'key_family'} };

    # Require that this actually have members!  For some reason, the grep just above is not blocking that problem...
    if (@key_families) {
        foreach my $key_family (@key_families) {
            my @gene_members = sort keys %{ $data_ref->{'family'}->{$key_family}->{'gene'} };
            my $gene_count   = @gene_members;
            my $annot_sum    = 0;

            # For each gene family, record its gene count:
            $data_ref->{'family'}->{$key_family}->{'gene_count'} = $gene_count;

            foreach my $gene_member (@gene_members) {
                my $annot = 0;
                if ( exists $data_ref->{'gene'}->{$gene_member}->{'annot'} ) {
                    $annot = $data_ref->{'gene'}->{$gene_member}->{'annot'};
                }
                else {
                    warn "Cannot map gene $gene_member to annot. count\n";
                }
                $annot_sum = $annot_sum + $annot;
            }
            $data_ref->{'family'}->{$key_family}->{'gene_count'} = $gene_count;
            $data_ref->{'family'}->{$key_family}->{'annot_sum'}  = $annot_sum;
        }

        # For each key gene, sort its families by annotation sum.
        @key_families = sort { $data_ref->{'family'}->{$a}->{'annot_sum'} <=> $data_ref->{'family'}->{$b}->{'annot_sum'} } @key_families;

        # For all the families that contain this gene, pick the one family that has the *lowest* aggregate annotation value.
        my $opt_fam = $key_families[0];
        $data_ref->{'final_gene'}->{$key_gene}->{'opt_fam'} = $opt_fam;

        # Get annot. sum, gene count, annot/gene value.
        my $opt_annot_sum = 0;
        if ( exists $data_ref->{'family'}->{$opt_fam}->{'annot_sum'} ) {
            $opt_annot_sum = $data_ref->{'family'}->{$opt_fam}->{'annot_sum'};
        }
        my $opt_gene_count = 0;
        if ( exists $data_ref->{'family'}->{$opt_fam}->{'gene_count'} ) {
            $opt_gene_count = $data_ref->{'family'}->{$opt_fam}->{'gene_count'};
        }
        my $anngene_ratio = 'n/a';
        if ( $opt_gene_count > 0 ) {
            $anngene_ratio = ($opt_annot_sum / $opt_gene_count);
        }

        # Link annot. sums to *genes*, to make it easier to sort them well later.
        $data_ref->{'final_gene'}->{$key_gene}->{'annot_sum'} = $opt_annot_sum;

        # Add annot/gene ratios to family annotations.
        $data_ref->{'family'}->{$opt_fam}->{'anngene_ratio'} = $anngene_ratio;
    }
}

# Ensure that @key_genes is nonredundant and universally has annot_sum values, then sort it by descending annotation levels.
@key_genes = grep { ( exists $data_ref->{'final_gene'}->{$_}->{'annot_sum'} ) } uniq(@key_genes);
@key_genes = sort { $data_ref->{'final_gene'}->{$a}->{'annot_sum'} <=> $data_ref->{'final_gene'}->{$b}->{'annot_sum'} } @key_genes;

# Finally, for each final key gene, print the relevant data.
foreach my $final_gene (@key_genes) {
    # This 'if' condition is crucial for weeding out genes that lack any opt_fam (which can happen if other taxa are being required).
    # Given that there is a defined $opt_fam, the other values ($gene_count, etc.) can be relied upon to exist.
    if ( exists $data_ref->{'final_gene'}->{$final_gene}->{'opt_fam'} ) {
        my $opt_fam       = $data_ref->{'final_gene'}->{$final_gene}->{'opt_fam'};
        my $gene_count    = $data_ref->{'family'}->{$opt_fam}->{'gene_count'};
        my $annot_sum     = $data_ref->{'family'}->{$opt_fam}->{'annot_sum'};
        my $anngene_ratio = $data_ref->{'family'}->{$opt_fam}->{'anngene_ratio'};

        # For $opt_fam, get a full membership list:
        my @fam_genes     = sort keys %{ $data_ref->{'family'}->{$opt_fam}->{'gene'} };
        my $fam_gene_list = join "; ", @fam_genes;

        # For $opt_fam, get a full genes-per-taxon count:
        my %tax_genes = ();
        foreach my $fam_gene (@fam_genes) {
            if ( $fam_gene =~ /\A ([a-z]+)\|\S+ \z/xms ) {
                my $taxon = $1;
                $tax_genes{$taxon}++;
            }
            else {
                die "Cannot parse taxon for gene $fam_gene in family $opt_fam\n";
            }
        }
        my @fam_taxa           = sort keys %tax_genes;
        my @fam_taxgene_counts = ();
        foreach my $fam_taxon (@fam_taxa) {
            # Avoid the barbarism of "(1 genes)".
            my $gene_word =	q{};
            if ( $tax_genes{$fam_taxon} >= 2 ) {
                $gene_word = 'genes';
            }
            else {
                $gene_word = 'gene';
            }
            my $fam_taxgene_count = "$fam_taxon ($tax_genes{$fam_taxon} $gene_word)";
            push @fam_taxgene_counts, $fam_taxgene_count;
        }
        my $fam_taxgene_text = join "; ", @fam_taxgene_counts;

        # Print header only once, at the start:
        print "$header\n" if $header;
        $header = q{};

        print "$final_gene";
        print "\t$opt_fam";
        print "\t$fam_taxgene_text";
        print "\t$fam_gene_list";
        print "\t$gene_count genes";
        print "\t$annot_sum";
        print "\t$anngene_ratio";
        print "\n";
    }
}

