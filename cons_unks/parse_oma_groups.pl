#!/usr/bin/env perl

# parse_oma_groups.pl -- Erich Schwarz <ems394@cornell.edu>, 10/5/2014.
# Purpose: given http://omabrowser.org/All/oma-groups.txt.gz and http://omabrowser.org/All/oma-uniprot.txt.gz, parse oma-groups.txt into OmaGroup_0001-to-UniProt (etc.).

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $uniprot_map = q{};
my $oma_groups  = q{};

my $help;

GetOptions ( 'uniprot_map=s' => \$uniprot_map,
             'oma_groups=s'  => \$oma_groups,
             'help'          => \$help,
);

if ( $help or (! $oma_groups ) or (! $uniprot_map ) ) {
    die "Format: parse_oma_groups.pl\n",
        "    --uniprot_map|-u  [oma-uniprot.txt table from OMA, linking OMA to UniProt protein IDs]\n",
        "    --oma_groups|-o   [oma-groups.txt table from OMA, with OMA IDs for proteins]\n",
        "    --help|-h         [print this message]\n",
        ;
}

open my $UNIPROT, '<', $uniprot_map;
while ( my $input = <$UNIPROT> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $oma_id  = $1;
        my $uniprot = $2;

        # OMA maps to both UniProt IDs and UniProt ACCs; we want ACCs only; ACCs lack underscores.  So...
        if ( $uniprot !~ /[_]/xms ) { 
            # We have to allow an OMA ID to map to multiple UniProt IDs, then sort out the chaos later.
            $data_ref->{'oma_id'}->{$oma_id}->{'uniprot_id'}->{$uniprot} = 1;
        }
    }
    elsif ( $input !~ /\A [#] /xms ) {
        die "From OMA/UniProt map table $uniprot_map, can't parse: $input\n";
    }
}
close $UNIPROT;

# This file needs to be opened twice.  The first time around, we are doing this *solely* so that I can pad group names with zeros in a well-controlled way later on.
open my $OMA, '<', $oma_groups;
while ( my $input = <$OMA> ) {
chomp $input;
    if ( $input =~ /\A (\d+) \t \S+ \t (.+) \z/xms ) {
        my $oma_group_no   = $1;
        if ( (! exists $data_ref->{'max_oma_group_no'} ) or ( $oma_group_no > $data_ref->{'max_oma_group_no'} ) ) { 
            $data_ref->{'max_oma_group_no' } = $oma_group_no;
        }
    }
}
close $OMA;

$data_ref->{'len_max_oma_group_no'} = length $data_ref->{'max_oma_group_no'};

open $OMA, '<', $oma_groups;
while ( my $input = <$OMA> ) { 

    # Sample imput lines:
    # [start of file]
    # # Orthologous groups from OMA release of Feb 2014
    # # This release has 862898 groups covering 6571118 proteins from 1613 species
    # # Format: group number<tab>Fingerprint<tab>tab-separated list of OMA Entry IDs
    # 1	GCQLRRS	ARATH02026	ARALY09035
    # 2	KDEHCES	ARATH05467	ARALY11434
    # 3	n/a	ECOLI04190	ECOD103527
    # [...]
    # 8	n/a	ECOLI01939	ECODH01875	ECOBW01756	ECOD101605
    # [...]
    # 23	ITRGNFF	CANAW01357	DEBHA03314	PICST00069	LODEL05715	DEKBR01151	CANTE03980	SPAPN01155
    # [...]
    # 862898	FANQKSH	WOLCO05133	PHYBL02392
    # [end of file]

    chomp $input;
    if ( $input =~ /\A (\d+) \t (\S+) \t (.+) \z/xms ) { 
        my $oma_group_no    = $1;
        my $oma_fingerprint = $2;
        my $oma_group_list  = $3;
        my @oma_group_seqs  = split /\t/, $oma_group_list;

        # Pad the group number with zeros so that the maximum number has none.
        my $sf_format = '%0' . $data_ref->{'len_max_oma_group_no'} . 'u';
        $oma_group_no = sprintf($sf_format, $oma_group_no) or die "Can't zero-pad OMA group number $oma_group_no\n";

        # make the 'n/a' fingerprint grep-friendier
        if ( $oma_fingerprint eq 'n/a' ) {
            $oma_fingerprint = 'no_fingerprint';
        }

        my $oma_group_id   = 'OMA_' . $oma_group_no . q{|} . $oma_fingerprint;
        foreach my $oma_id (@oma_group_seqs) { 
            if ( exists $data_ref->{'oma_id'}->{$oma_id}->{'uniprot_id'} ) {
                my @uniprot_ids = sort keys %{ $data_ref->{'oma_id'}->{$oma_id}->{'uniprot_id'} };
                foreach my $uniprot_id (@uniprot_ids) {
                    print "$oma_group_id\t$uniprot_id\n";
                }
            }
        }
    }
}
close $OMA;

