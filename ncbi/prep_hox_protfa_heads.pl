#!/usr/bin/env perl

# prep_hox_protfa_heads.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/22/2008.
# Purpose: given homology texts and Hox paper sequence names, prepare protein FASTA headers for Sequin.

use strict;
use warnings;

my $input_argline = join ' ', @ARGV;
if ($input_argline !~ /\A --homologies \s \S.+\S \s --fastas \s \S.+\S \z /xms ) { 
    die "Format: ./prep_hoxfa_heads.pl --homologies [table(s)] --fastas [FASTA(s)]\n";
}
if ($input_argline =~ /\A --homologies \s (\S.+\S) \s --fastas \s (\S.+\S) \z /xms ) { 
    my $homology_list = $1;
    my $fasta_list    = $2;
    my %id2annot      = ();
    my @homolog_files = split /\s/, $homology_list;
    my @fasta_files   = split /\s/, $fasta_list;
    foreach my $hom_file (@homolog_files) { 
        if (! -e $hom_file) { 
            die "Can't find putative homology file $hom_file!\n";
        }
        open my $HOMS, '<', $hom_file or die "Can't open homology file $hom_file!\n";

# Typical input lines:
# Cbre_JD01.003   422 aa  WBGene00016653|C44E4.4 [*] (1 elegans, 1 briggsae, 1 remanei, 1 brenneri).
# Cbre_JD01.004   4,217 aa        WBGene00016650|C44E4.1 and WBGene00016656|C44E4.7 (2 elegans ...
# 
# Csp3_JD04.010   211 aa  WBGene00001174|egl-5 [*] (1 elegans, 1 briggsae, 1 remanei, 1 ps1010).
# Csp3_JD04.011   1,086 aa        WBGene00000768|cor-1 and WBGene00007983|C36E8.4 (2 elegans

        LOOP1:
        while ( my $input = <$HOMS> ) {
            chomp $input;
            if ( $input =~ / \A 
                             ( C\w+_JD\d+\.\d+ ) 
                             \s+ \S+ \s aa \s+ 
                             ( WBGene \d+ .+ ) 
                             \( \d+ \s elegans 
                           /xms ) { 
                my $prot_id1 = $1;
                my $annot1   = $2;
                $annot1 =~ s/(\A\s+|\s+\z)//g;
                if ($annot1 =~ / \A (WBGene\d+\|\S+) 
                                \s+ \[\*\] \z 
                             /xms ) { 
                    $annot1 = "[prot_desc=Orthologous to C. elegans $1.]";
                    $id2annot{$prot_id1} = $annot1;
                    next LOOP1;
                }
                $annot1 = "[prot_desc=Similar to C. elegans $annot1.]";
                $id2annot{$prot_id1} = $annot1;
            }
            elsif ( $input =~ / \A ( C\w+_JD\d+\.\d+ ) /xms ) { 
                my $prot_id2 = $1;
                my $annot2   = '[prot_desc=Hypothetical protein.]';
                $id2annot{$prot_id2} = $annot2;
            }
        }
        # End of LOOP1.
        close $HOMS or die "Can't close filehandle to $hom_file!\n";
    }
    foreach my $fa_file (@fasta_files) { 
        if (! -e $fa_file) {
            die "Can't find putative FASTA file $fa_file!\n";
        }
        open my $FASTA, '<', $fa_file or die "Can't open FASTA file $fa_file!\n";
        while ( my $input = <$FASTA> ) { 
            chomp $input;
            if ( $input =~ /\A > (([^\.]+)\.\d+) /xms ) { 
                my $fa_id      = $1;
                my $ncbi_fa_id = $2;

                # Insert NCBI's 'lcl|' boilerplate in front of contig and gene IDs;
                #    use NCBI's strange system of having the contig name after '>'.
                print ">lcl|$ncbi_fa_id";
                print "  [protein=$fa_id]  [gene=lcl|$fa_id]";
                print "  $id2annot{$fa_id}" if (exists $id2annot{$fa_id});

                print "\n";
            }
            else { 
                print "$input\n";
            }
        }
        close $FASTA or die "Can't close filehandle to $fa_file!\n";
    }
}

