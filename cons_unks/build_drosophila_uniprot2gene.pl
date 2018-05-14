#!/usr/bin/env perl

use strict;
use warnings;

my $data_ref;

my $id_table_1 = $ARGV[0];
my $db_name    = $ARGV[1];
my $id_table_2 = $ARGV[2];

open my $ID1, '<', $id_table_1 or die "Can't open MOD ID table 1 $id_table_1: $!";
while (my $input = <$ID1>) { 
    chomp $input;
    # Sample input:
    # P32234	FlyBase	FBgn0010339
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) { 
        my $uniprot = $1;
        my $db_ref  = $2;
        my $mod_id  = $3;
        if ( $db_ref eq $db_name ) {
            if ( $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} ) {
                my @known_mod_ids = sort keys %{$data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
                push @known_mod_ids, $mod_id;
                my $known_mod_id_text = join ' versus ', @known_mod_ids;
                warn "Multiple mapping of UniProt $uniprot to gene IDs $known_mod_id_text\n";
            }
            if ( $mod_id !~ /\A FBgn\d+ \z/xms ) { 
                die "Can't make sense of putative MOD ID $mod_id in input from MOD ID table 1 $id_table_1: $input\n";
            }
            # Note: we have to do this, because there are instances of *identical* proteins encoded by two or more genes...
            $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'}->{$mod_id} = 1;
        }
    }
    else { 
        die "Can't parse input from MOD ID table 1 $id_table_1: $input\n";
    }
}
close $ID1 or die "Can't close filehandle to MOD ID table 1 $id_table_1: $!";

open my $ID2, '<', $id_table_2 or die "Can't open MOD ID table 2 $id_table_2: $!";
while (my $input = <$ID2>) {
    chomp $input;
    if ( ( $input !~ /\A [#] /xms ) and ( $input =~ /\A  (\S+) \t (FBgn\d+) \t [^\t]* \t [^\t]* \t \S+ /xms ) ) { 
        # Sample input:
        # ## Using datasource: dbi:Pg:dbname=fb_2013_04_reporting;host=flysql9;port=5432...
        # 
        # ##gene_symbol	primary_FBgn#	nucleotide_accession	na_based_protein_accession	UniprotKB/Swiss-Prot/TrEMBL_accession
        # d	FBgn0262029	AC004728		
        # d	FBgn0262029	AE014134	ADV36979	
        # d	FBgn0262029			Q8T6L9

        my $common_name = $1;
        my $mod_id      = $2;
        # The last \S+ should be UniProt, though for simplicity we'll ignore it.

        my $long_name   = $mod_id . q{|} . $common_name;

        # Check this *first*, or it will always get triggered...
        if ( ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ) and ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ne $long_name ) ) {
            die "Inconsistent mapping of MOD ID $mod_id to long gene names $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} and $long_name\n";
        }

        # Finally:
        $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} = $long_name;
    }
    elsif ( ( $input !~ /\A # /xms ) and ( $input !~ /\A \s* \z/xms ) ) {
        die "Can't parse input from MOD ID table 2 $id_table_2: $input\n";
    }
}
close $ID2 or die "Can't close filehandle to MOD ID table 2 $id_table_2: $!";

my @uniprots = sort keys %{ $data_ref->{'uniprot'} };
my $uniprot_count = @uniprots;
foreach my $uniprot (@uniprots) {
    my @mod_ids = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'mod_gene_id'} };
    foreach my $mod_id (@mod_ids) {
        if ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ) { 
            my $long_name = $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'};
            print "$uniprot\t$long_name\n";
        }
        else {
            warn "Can't find long name for uniprot $uniprot and MOD ID $mod_id\n";
        }
    }
}

