#!/usr/bin/env perl

# map_various_orths_to_ref_density.pl -- Erich Schwarz <ems394@cornell.edu>, 10/4/2014.
# Purpose: use UniProt-[various orthology DB] mappings, with other tables, to get a ranking of conserved PantherDB groups by known-ness or un-known-ness.
# Note: For instance, the UniProt-[various orthology DB] mappings might, themselves, to be extracted from *.dat files, preferably organism-specific; or they might be from an orthology database such as TreeFam or OMA.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use autodie;

my $data_ref;

my @orthogroups2uniprots = ();
my @uniprots2orthogroups = ();
my @uniprots2genes       = ();
my @genes2densities      = ();

my $help;

GetOptions ( 'orthogroups2uniprots|o=s{,}' => \@orthogroups2uniprots,
             'uniprots2orthogroups|u=s{,}' => \@uniprots2orthogroups,
             'uniprots2genes|g=s{,}'       => \@uniprots2genes,
             'genes2densities|d=s{,}'      => \@genes2densities,
             'help'                      => \$help, 
);


if ( $help or ( (! @orthogroups2uniprots) and (! @uniprots2orthogroups ) ) or (! @uniprots2genes) or (! @genes2densities) ) { 
    die "Format: map_various_orths_to_ref_density.pl\n",
        "    --orthogroups2uniprots|-o   [1+ precomputed [orthogroup]-to-UniProt tables]\n",
        "    --uniprots2orthogroups|-u   [1+ precomputed UniProt-to-[orthogroup] tables]\n",
        "    [There can be mixed sets of -o and -u table files, but there *must* be at least one table file of at least one table type.]\n",
        "    --uniprots2genes|-g         [1+ uniprot to gene tables; if using several, make sure no ambiguous short names!]\n",
        "    --genes2densities|-d        [1+ gene to annotation density tables; again, avoid short ambiguous names]\n",
        "    --help|-h                   [print this message]\n",
        ;
}

foreach my $orthogroup2uniprot (@orthogroups2uniprots) { 
    open my $OGROUP2UPROT, '<', $orthogroup2uniprot;
    while (my $input = <$OGROUP2UPROT>) {
        chomp $input;

        # Sample input:
        # eggNOG|COG5040  P31946
        # InParanoid|P31946       P31946
        # OMA|MGREYRE     P31946
        # [etc.]

        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $orthogroup = $1;
            my $uniprot    = $2;
            $data_ref->{'orthogroup'}->{$orthogroup}->{'uniprot'}->{$uniprot} = 1;
         }
         else { 
             # Allow comments to be passed over in silence if they are marked with '# ' at their beginning;
             #    otherwise, die loudly.
             if (  $input !~ /\A [#] /xms ) { 
                 die "From orthogroup2uniprot file $orthogroup2uniprot, can't parse: $input\n";
             }
         }
    }
    close $OGROUP2UPROT;
}

foreach my $uniprot2orthogroup (@uniprots2orthogroups) {
    open my $UPROT2OGROUP, '<', $uniprot2orthogroup;
    while (my $input = <$UPROT2OGROUP>) {
        chomp $input;
             
        # Sample input:
        # eggNOG|COG5040  P31946
        # InParanoid|P31946       P31946
        # OMA|MGREYRE     P31946
        # [etc.]
    
        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $uniprot    = $1;
            my $orthogroup = $2;
            $data_ref->{'orthogroup'}->{$orthogroup}->{'uniprot'}->{$uniprot} = 1;
         }
         else {
             # Allow comments to be passed over in silence if they are marked with '# ' at their beginning;
             #    otherwise, die loudly.
             if (  $input !~ /\A [#] /xms ) {
                 die "From uniprot2orthogroup file $uniprot2orthogroup, can't parse: $input\n";
             }
         }
    }
    close $UPROT2OGROUP;
}

foreach my $uniprot2gene (@uniprots2genes) {
    open my $GENE, '<', $uniprot2gene;
    while (my $input = <$GENE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $uniprot = $1;
            my $gene    = $2;
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
        }
    }
    close $GENE;
}

foreach my $gene2density (@genes2densities) { 
    open my $DENSITY, '<', $gene2density;
    while (my $input = <$DENSITY>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $gene    = $1;
            my $density = $2;
            $data_ref->{'gene'}->{$gene}->{'density'} = $density;
        }
    }
    close $DENSITY;
}

my @orthogroups = sort keys %{ $data_ref->{'orthogroup'} };
my @orthogroup_reports = ();

foreach my $orthogroup (@orthogroups) { 
    my $total_ref_density = 0;
    my @genes             = ();
    my @uniprots = sort keys %{ $data_ref->{'orthogroup'}->{$orthogroup}->{'uniprot'} };
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
    }

    # No more than four decimal places, please!  But put this rounding *after* all the effects of all the genes have been summed up.
    $total_ref_density = sprintf "%.4f", $total_ref_density;

    push @orthogroup_reports , "$orthogroup\t$total_ref_density\t$gene_count\t$gene_list";
}

@orthogroup_reports = sort { rank_by_annot($b) <=> rank_by_annot($a) } @orthogroup_reports;

foreach my $orthogroup_report (@orthogroup_reports) { 
    print "$orthogroup_report\n";
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

