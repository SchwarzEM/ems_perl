#!/usr/bin/env perl

# stitch_uniprot_amiD_data_28jul2014.pl -- Erich Schwarz, 7/28/2014.
# Purpose: link data from indiv_proteome_files/uniprot_sprot.07jul2014.fa or indiv_proteome_files/uniprot_trembl.07jul2014.fa to pp_secondary_amiD_domain_subset_25jul2014a.fa.

use strict;
use warnings;
use autodie;

my @uniprot_seqs   = @ARGV;
my $prefixed_amids = pop @uniprot_seqs;

my %pref2class = (
    Arch_ => 'Archaea',
    Arth_ => 'Arthropods',
    Bact_ => 'Bacteria',
    Euka_ => 'Eukaryotes',
    Metg_ => 'Metagenome',
    Mzoa_ => 'Metazoa',
    Nema_ => 'Nematodes',
    Vert_ => 'Vertebrates',
    Viru_ => 'Viruses',
);

my $data_ref;

open my $PREF, '<', $prefixed_amids;
while (my $input = <$PREF>) {
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ( ([A-Z][a-z]+ [_]) ( [A-Za-z0-9]+ [.]* [_] [^\s\/]+ ) ) \/ \d+ [-] \d+ /xms ) { 
            my $two_prefixed_name = $1;
            my $class_prefix      = $2;
            my $one_prefixed_name = $3;
            if (     ( exists $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'}                ) 
                 and ( $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'} ne $two_prefixed_name ) ) {
                my $new_two_pr_name = $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'};
                die "One-prefixed name $one_prefixed_name has two different two-prefixed names ($two_prefixed_name vs. $new_two_pr_name) -- see header: $input\n";
            }
            $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'} = $two_prefixed_name;
            my $class = $class_prefix;
            $class    = $pref2class{$class};
            $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'class'} = $class;
        }
        else { 
            die "From prefixed amiD file $prefixed_amids, can't parse header: $input\n";
        }
    }
}
close $PREF;

foreach my $uniprot_seq (@uniprot_seqs) {
    open my $UNIPROT, '<', $uniprot_seq;
    while (my $input = <$UNIPROT>) {
        chomp $input;
        if ( $input =~ /\A > /xms ) {
            if ( $input =~ /\A > ( [^_\s]+ [_] (\S+) ) \t [^\t]* \t [^\t]* \t ([^\t]*) \t ([^\t]*) \t NCBI_TaxID=\d+ \t uniprot [_] [a-z]+ \.dat /xms ) { 
                my $one_prefixed_name = $1;
                my $orig_name         = $2;
                my $gene              = $3;
                my $species           = $4;
                if (     ( exists $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'class'}             ) 
                     and ( exists $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'} ) ) { 
                    my $class = $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'class'};
                    my $two_prefixed_name = $data_ref->{'one_prefixed_name'}->{$one_prefixed_name}->{'two_prefixed_name'};
                    print "$two_prefixed_name\t$class\t$species\t$orig_name\t$gene\n";
                }
            }
            else {
                die "From UniProt file $uniprot_seq, can't parse header: $input\n";
            }
        }
    }
    close $UNIPROT;
}

