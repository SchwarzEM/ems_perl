#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $tx2gene2spp = q{};
my $ethresh     = 1;
my @blastps     = ();

my $data_ref;

my $help;

GetOptions ( 'blastps=s{,}'  => \@blastps,
             'tx2gene2spp=s' => \$tx2gene2spp,
             'ethresh=s'     => \$ethresh,
             'help'          => \$help,   );

if ( (! looks_like_number($ethresh) ) or ( $ethresh < 0 ) or ( $ethresh > 1 ) ) {
    warn "E-value threshold is not a number between 0 and 1: $ethresh\n";
}

if ($help or (! @blastps ) or (! $tx2gene2spp ) or (! looks_like_number($ethresh) ) ) {
    die "Format: parse_orphan_blastps_29nov2016.pl\n",
        "    --blastps|-b      [BlastP files]\n",
        "    --tx2gene2spp|-t  [transcript to gene to species table]\n",
        "    --ethresh|-e      [E value threshold for significant hits, between 0 and 1]\n",
        "    --help|-h         [print this message]\n",
        ;
}

open my $TX2GENE2SPP, '<', $tx2gene2spp;
while (my $input = <$TX2GENE2SPP>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
        my $prot    = $1;
        my $gene    = $2;
        my $species = $3;

        $data_ref->{'prot'}->{$prot}->{'gene'}    = $gene;
        $data_ref->{'prot'}->{$prot}->{'species'} = $species;
    }
    else {
        die "In tx2gene2spp table ($tx2gene2spp), cannot parse: $input\n";
    }
}
close $TX2GENE2SPP;

foreach my $blastp (@blastps) {
    open my $BLASTP, '<', $blastp;
    while (my $input = <$BLASTP>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) {
            my $prot  = $1;
            my $match = $2;
            my $eval  = $3;

            if (! exists $data_ref->{'prot'}->{$prot}->{'species'} ) {
                die "In BlastP report ($blastp), cannot find species for query protein $prot\n";
            }

            if (! exists $data_ref->{'prot'}->{$match}->{'species'} ) {
                die "In BlastP report ($blastp), cannot find species for match protein $match\n";
            }

            if (! looks_like_number($eval) ) {
                die "In BlastP report ($blastp), non-numerical e-value in: $input\n";
            }

            my $query_species = $data_ref->{'prot'}->{$prot}->{'species'};
            if ( $query_species ne 'nigoni' ) {
                die "In BlastP report ($blastp), query protein $prot is non-nigoni species $query_species!\n";
            }

            my $match_species = $data_ref->{'prot'}->{$match}->{'species'};

            if ( ( $eval <= $ethresh ) and ( $match_species ne 'nigoni' ) ) {
                my $gene = $data_ref->{'prot'}->{$prot}->{'gene'};
                $data_ref->{'non_orphan'}->{$gene}->{'other_species'}->{$match_species} = 1;
            }
        }
        else {
            die "In BlastP report ($blastp), cannot parse: $input\n";
        }
    }
    close $BLASTP;
}

my @non_orphan_genes = sort keys %{ $data_ref->{'non_orphan'} };

foreach my $non_orphan_gene (@non_orphan_genes) {
    my @non_orphan_species      = sort keys %{ $data_ref->{'non_orphan'}->{$non_orphan_gene}->{'other_species'} };
    my $non_orphan_species_text = join '; ', @non_orphan_species;
    print "$non_orphan_gene\t$non_orphan_species_text\n";
}
