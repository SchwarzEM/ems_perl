#!/usr/bin/env perl

# map_panther_to_ref_density.pl -- Erich Schwarz <ems394@cornell.edu>, originally 8/26/2013; small but significant revision on 10/3/2014, to avoid overrounding.
# Purpose: use precomputed UniProt-PantherDB mappings, with other tables, to get a ranking of conserved PantherDB groups by known-ness or un-known-ness.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $data_ref;

my @panthers2uniprots = ();
my @uniprots2genes  = ();
my @genes2densities = ();

my $help;

GetOptions ( 'panthers2uniprots=s{,}' => \@panthers2uniprots,
             'uniprots2genes=s{,}'    => \@uniprots2genes,
             'genes2densities=s{,}'   => \@genes2densities,
             'help'                   => \$help, 
);


if ( $help or (! @panthers2uniprots) or (! @uniprots2genes) or (! @genes2densities) ) { 
    die "Format: map_panther_to_ref_density.pl\n",
        "    --panthers2uniprots|-p   [1+ precomputed UniProt/PantherDB tables from PantherDB.org]\n",
        "    --uniprots2genes|-u      [1+ uniprot to gene tables; if using several, make sure no ambiguous short names!]\n",
        "    --genes2densities|-g     [1+ gene to annotation density tables; again, avoid short ambiguous names]\n",
        "    --help|-h                [print this message]\n",
        ;
}

foreach my $panther2uniprot (@panthers2uniprots) { 
    open my $PANTHER, '<', $panther2uniprot or die "Can't open panther2uniprot file $panther2uniprot: $!";
    while (my $input = <$PANTHER>) {
        chomp $input;
        # Sample input:
        # CAEEL|WB=WBGene00008634|UniProtKB=Q8I4L6		PTHR21447:SF0	UNCHARACTERIZED	SUBFAMILY NOT NAMED
        # DROME|FB=FBgn0028375|UniProtKB=Q9V3R8		PTHR13929:SF0	1,4-DIHYDROXY-2-NAPHTHOATE OCTAPRENYLTRANSFERASE	SUBFAMILY NOT NAMED	transferase activity#GO:0016740;catalytic activity#GO:0003824	lipid metabolic process#GO:0006629;vitamin biosynthetic process#GO:0009110;vitamin metabolic process#GO:0006766;primary metabolic process#GO:0044238;metabolic process#GO:0008152		transferase#PC00220	
        if ( $input !~ /\A [#] /xms ) { 
            if ( $input =~ /\A \S+ \| UniProtKB = (\S+) \t [^\t]* \t (PTHR\d+) : SF\d+ \t ([^\t]+) \t /xms ) { 
                my $uniprot      = $1;
                my $panther_id   = $2;
                my $panther_text = $3;

                # Have to deal with PantherDB free-text, which includes ' ' and ',' at least.
                $panther_text =~ s/\s/_/g;
                $panther_text =~ s/,/./g;

                my $panther_name = $panther_id . q{|} . $panther_text;
                $data_ref->{'panther'}->{$panther_name}->{'uniprot'}->{$uniprot} = 1;
            }
            else { 
                die "From panther2uniprot file $panther2uniprot, can't parse: $input\n";
            }
        }
    }
    close $PANTHER or die "Can't close filehandle to panther2uniprot file $panther2uniprot: $!";
}

foreach my $uniprot2gene (@uniprots2genes) {
    open my $GENE, '<', $uniprot2gene or die "Can't open uniprot2gene file $uniprot2gene: $!";
    while (my $input = <$GENE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $uniprot = $1;
            my $gene    = $2;
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
        }
    }
    close $GENE or die "Can't close filehandle to uniprot2gene file $uniprot2gene: $!";
}

foreach my $gene2density (@genes2densities) { 
    open my $DENSITY, '<', $gene2density or die "Can't open gene2density file $gene2density: $!";
    while (my $input = <$DENSITY>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $gene    = $1;
            my $density = $2;
            $data_ref->{'gene'}->{$gene}->{'density'} = $density;
        }
    }
    close $DENSITY or die "Can't close filehandle to gene2density file $gene2density: $!";
}

my @panther_domains = sort keys %{ $data_ref->{'panther'} };
my @panther_reports = ();
foreach my $panther (@panther_domains) { 
    my $total_ref_density = 0;
    my @genes             = ();
    my @uniprots = sort keys %{ $data_ref->{'panther'}->{$panther}->{'uniprot'} };
    foreach my $uniprot (@uniprots) {
        my @new_genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} };
        push @genes, @new_genes;
    }

    @genes         = sort @genes;
    @genes         = uniq @genes;
    my $gene_list  = join '; ', @genes;
    my $gene_no    = @genes;
    my $gene_count = q{};
    if ( $gene_no >= 2 ) { 
        $gene_count = "$gene_no genes";
    }
    if ( $gene_no == 1 ) {
        $gene_count = "$gene_no gene";
    }

    foreach my $gene (@genes) { 
        my $gene_ref_weight = 0;

        if ( exists $data_ref->{'gene'}->{$gene}->{'density'} ) { 
            $gene_ref_weight = $data_ref->{'gene'}->{$gene}->{'density'};
        }

        $total_ref_density = $total_ref_density + $gene_ref_weight;

        # Rounding was done here -- repeatedly -- in the Aug. 2013 version of this code.
        # In retrospect, that seems like a great way to have imprecision; and in any event, it is needless.
    }

    # No more than four decimal places, please!
    # However, do the rounding only *once*, *after* all the gene annotation densities have been summed for the homology group.
    $total_ref_density = sprintf "%.4f", $total_ref_density;

    push @panther_reports, "$panther\t$total_ref_density\t$gene_count\t$gene_list";
}

@panther_reports = sort { rank_by_annot($b) <=> rank_by_annot($a) } @panther_reports;

foreach my $panther_report (@panther_reports) { 
    print "$panther_report\n";
}

sub rank_by_annot {
    my $_report = $_[0];
    if ( $_report =~ /\A \S+ \t (\S+) \b/xms ) { 
        my $_annot_density = $1;
        return $_annot_density;
    }
    else { 
        die "Subroutine rank_by_annot can't parse annotation report: $_report\n";
    }
}

