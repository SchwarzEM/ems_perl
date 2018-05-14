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
    # Q5U3F8  ZFIN    ZDB-GENE-041114-89
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
            if ( $mod_id !~ /\A ZDB[-]GENE[-]\d+[-]\d+ \z/xms ) { 
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
    # if ( $input =~ /\A (\S+) \t SO:0000704 \t (\S+) /xms ) {
    if ( $input =~ /\A (ZDB[-]GENE[-]\d+[-]\d+) \t SO:0000704 \t ([^\t]*) \t \S+ \t /xms ) { 
        # Sample input:
        # ZDB-GENE-000112-32      SO:0000704      couptf3 O42533  
        # ZDB-GENE-000112-34      SO:0000704      couptf4 O42534  
        # ZDB-GENE-000112-36      SO:0000704      couptf5 O42535  
        my $mod_id      = $1;
        my $common_name = $2;
        my $long_name   = $mod_id ;
        if ( ( $common_name =~ /\S/xms ) and ( $common_name ne $mod_id ) ) {
            $long_name = $long_name . q{|} . $common_name;
        }

        # Check this *first*, or it will always get triggered...
        if ( ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ) and ( $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} ne $long_name ) ) {
            die "Inconsistent mapping of MOD ID $mod_id to long gene names $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} and $long_name\n";
        }

        # Finally:
        $data_ref->{'mod_gene_id'}->{$mod_id}->{'long_name'} = $long_name;
    }
    else {
        die "Must learn to parse: $input\n";
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

