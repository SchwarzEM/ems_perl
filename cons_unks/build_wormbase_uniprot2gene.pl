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
    # P41932	WormBase	WBGene00003920
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
            if ( $mod_id !~ /\A WBGene\d+ \z/xms ) { 
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
    if ( $input =~ /\A > (\S+) \s+ CE\d+ \s+ (WBGene\d+) \s+ (\S+) /xms ) { 
        # Sample input:
        # >2L52.1	CE32090	WBGene00007063	Zinc finger, C2H2 type	status:Partially_confirmed	UniProt:A4F336	protein_id:CCD61130.1
        # >3R5.1a	CE24758	WBGene00007065	locus:pot-3	status:Partially_confirmed	UniProt:G5EFG7	protein_id:CAA21777.2

        my $cds_name    = $1;
        my $mod_id      = $2;
        my $common_name = $3;
        $cds_name       =~ s/[a-z]+\z//;
        my $long_name   = $mod_id . q{|} . $cds_name;

        if ( $common_name !~ /\A locus : /xms ) { 
            $common_name = q{};
        }
        if ( $common_name =~ /\A locus : (\S+) /xms ) {
            $common_name = $1;
            $long_name   = $long_name . q{|} . $common_name;
        }

        # Check this *first*, or it will always get triggered...
        if ( ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ) and ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ne $long_name ) ) {
            die "Inconsistent mapping of MOD ID $mod_id to long gene names $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} and $long_name\n";
        }

        # Finally:
        $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} = $long_name;
    }
    elsif ( $input =~ /\A > /xms ) {
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

