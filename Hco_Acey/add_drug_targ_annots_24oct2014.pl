#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);

my $acey_pfams     = $ARGV[0];
my $acey_interpros = $ARGV[1];
my $cel_pfams      = $ARGV[2];
my $cel_interpros  = $ARGV[3];
my $drug_annots    = $ARGV[4];
my $big_data_table = $ARGV[5];

my $data_ref;

open my $ACEY_PFAM, '<', $acey_pfams;
while (my $input = <$ACEY_PFAM>) { 
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t [^\t]+ \t (Acey_\S+) \z/xms ) { 
        my $motif = $1;
        my $gene  = $2;
        # Correct retro-nomenclature to current nomenclature.
        $gene =~ s/Acey_2012.08.05_/Acey_s/;
        $data_ref->{'motif'}->{$motif}->{'acey_gene'}->{$gene} = 1;
    }
    else {
        if ( $input !~ /\A Accession [ ] no\. /xms ) { 
            die "From Acey PFAMs file $acey_pfams, cannot parse: $input\n";
        }
    }
}
close $ACEY_PFAM;

open my $ACEY_IPRO, '<', $acey_interpros;
while (my $input = <$ACEY_IPRO>) {
    chomp $input;
    # Sample input line:
    # HMMPanther	PTHR10169:SF19	DNA TOPOISOMERASE 2	Acey_2012.08.05_0064.g3495
    if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t [^\t]+ \t (Acey_\S+) \z/xms ) {
        my $motif_cat = $1;
        my $motif_acc = $2;
        my $gene      = $3;
        my $motif     = $motif_cat . q{|} . $motif_acc;
        # Correct retro-nomenclature to current nomenclature.
        $gene =~ s/Acey_2012.08.05_/Acey_s/;
        $data_ref->{'motif'}->{$motif}->{'acey_gene'}->{$gene} = 1;
    }
    else {
        if ( $input !~ /\A Subsidiary [ ] database \t Accession[ ] no\. /xms ) {
            die "From Acey InterPro file $acey_interpros, cannot parse: $input\n";
        }
    }
}
close $ACEY_IPRO;

open my $CEL_PFAM, '<', $cel_pfams;
while (my $input = <$CEL_PFAM>) {
    chomp $input;
    # Sample input line:
    # PF01274.17 Malate_synthase WBGene00001564|C05E4.9|icl-1
    if ( $input =~ /\A ([^\t]+) \t [^\t]+ \t WBGene\d+ \S* \| (\S+?) \z/xms ) {
        my $motif = $1;
        my $gene  = $2;
        $data_ref->{'cel_gene'}->{$gene}->{'motif'}->{$motif} = 1;
    }
    else {
        if ( $input !~ /\A PFAM [ ] accession [ ] no\. /xms ) {
            die "From Cel PFAMs file $cel_pfams, cannot parse: $input\n";
        }
    }
}
close $CEL_PFAM;

open my $CEL_IPRO, '<', $cel_interpros;
while (my $input = <$CEL_IPRO>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t [^\t]+ \t WBGene\d+ \S* \| (\S+?) \z/xms ) {
        my $motif_cat = $1;
        my $motif_acc = $2;
        my $gene      = $3;
        my $motif     = $motif_cat . q{|} . $motif_acc;
        $data_ref->{'cel_gene'}->{$gene}->{'motif'}->{$motif} = 1;
    }
    else {
        if ( $input !~ /\A Subsidiary [ ] database \t Accession[ ] no\. /xms ) {
            die "From Cel InterPro file $cel_interpros, cannot parse: $input\n";
        }
    }
}
close $CEL_IPRO;

open my $DRUG_ANNOTS, '<', $drug_annots;
while (my $input = <$DRUG_ANNOTS>) {
    chomp $input;

    # Sample input lines:
    # 4-coumarate:coenzyme A ligase, class I	10	acs-10	n/a	pmid10417722
    if ( $input =~ /\A ([^\t]+) \t \d+ \t ([^\t]+) \t ([^\t]+) \t [^\t]+ \z/xms ) { 
        my $target_type   = $1;
        my $cel_gene_list = $2;
        my $drug_type     = $3;
 
        if ( $drug_type ne 'n/a' ) {
            $drug_type = "; $drug_type";
        }
        if ( $drug_type eq 'n/a' ) { 
            $drug_type = q{};
        }
        my $drug_annot = "\"$target_type$drug_type\"";

        my @cel_genes = split /,\s+/, $cel_gene_list;
        my @acey_genes = ();
        foreach my $cel_gene (@cel_genes) { 
            my @mots = sort keys %{ $data_ref->{'cel_gene'}->{$cel_gene}->{'motif'} };
            foreach my $motif (@mots) { 
                my @new_acey_genes = sort keys %{ $data_ref->{'motif'}->{$motif}->{'acey_gene'} };
                push @acey_genes, @new_acey_genes;
            }
        }
        @acey_genes = sort @acey_genes;
        @acey_genes = uniq @acey_genes; 
        foreach my $acey_gene (@acey_genes) {
            if ( exists $data_ref->{'acey_gene'}->{$acey_gene}->{'drug_annot'} ) { 
                die "For Acey gene $acey_gene, two annotations: $drug_annot and $data_ref->{'acey_gene'}->{$acey_gene}->{'drug_annot'}\n";
            }
            else {
                $data_ref->{'acey_gene'}->{$acey_gene}->{'drug_annot'} = $drug_annot;
            }
        }
    }
    else { 
        if ( $input !~ /\A Protein [ ] class \t /xms ) { 
            die "From drug annots file $drug_annots, cannot parse: $input\n";
        }
    }
}
close $DRUG_ANNOTS;

open my $BIG_DATA, '<', $big_data_table;
while (my $input = <$BIG_DATA>) {
    chomp $input;
    if ( $input =~ /\A ([^\t]+) \t (?: [^\t]* \t){25} ([^\t]*) \t /xms ) { 
        my $gene             = $1;
        my $drug_status      = $2;
        my $prev_drug_status = q{};
        if ( $gene =~ /\A Acey_s\d+\.g\d+ \z/xms ) { 
            # Categorically reject *any* input line with two (or more) instances of the plaintext string 'Drug_target':
            if ( $input =~ / Drug_target .* Drug_target /xms ) {
                die "Two occurrences of text \"Drug_target\" in data line; $input\n";
            }

            # Edit out any unsupported annotations in the big data table; record *all* incongruities:
            elsif ( ( ( exists $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'} ) and ( $drug_status ne 'Drug_target' ) ) 
                 or
                 ( (! exists $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'} ) and ( $drug_status eq 'Drug_target' ) ) 
               ) {

                # Record the incongruity as an error:
                warn "Incongruity between annotation and big data table for: $input\n";
                warn "Big data drug status is: \"$drug_status\"\n";
                if ( exists $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'} ) {
                    $prev_drug_status = $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'};
                }
                warn "Recorded annotation is: \"$prev_drug_status\"\n";

                # Edit out the unsupported annotation:
                $input =~ s/Drug_target//;
            }

            # Enhance any nonpathological drug target annotations found in the big data file:
            elsif ( ( exists $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'} ) and ( $drug_status eq 'Drug_target' ) ) {
                if ( $input =~ / Drug_target .* Drug_target /xms ) { 
                    die "Two occurrences of text \"Drug_target\" in data line; $input\n";
                }
                my $drug_annot = 'Drug_target|' . $data_ref->{'acey_gene'}->{$gene}->{'drug_annot'};
                $input =~ s/Drug_target/$drug_annot/;
            }

            # Enforce successful parsing of any line that has a non-zero $drug_status value:
            elsif ($drug_status) {
                die "From big data table $big_data_table, cannot parse: $input\n";
            }
        }
    }
    print "$input\n";
}
close $BIG_DATA;

