#!/usr/bin/env perl

# genes2omcls.pl -- Erich Schwarz <ems394@cornell.edu>, 2/13/2013; generalized to deal with OrthoMCL-ized OrthoFinder data on 9/13/2016; updated headers on 7/29/2020.
# Purpose: given OrthoMCL-ized OrthoFinder (or real OrthoMCL) report text (file or piped) and a named taxon, annotate each taxon's genes/proteins with its OrthoMCL groups; optionally, summarize OrthoMCL.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles       = ();
my $taxon         = q{};
my $maximum       = 0;
my $summary;
my $help;

my $data_ref;

GetOptions ( 'infile=s{,}' => \@infiles,
             'taxon=s'     => \$taxon,
             'summary'     => \$summary,
             'max=i'       => \$maximum,
             'help'        => \$help,    );

if ( $help or (! $taxon ) ) {
    die "Format: genes2omcls.pl\n",
        "        --summary|-s  [summarize the OrthoMCL groups]\n",
        "        --taxon|-t    [taxon to count]\n", 
        "        --infile|-i   [input files or '-' for stream]\n",
        "        --max|-m      [set max. characters in ortholog description text; e.g., something less than 32,767 is useful for not killing Excel imports]\n",
        "        --help|-h     [print this message]\n",
        ;
}

if ( (! looks_like_number($maximum) ) or ( $maximum != int($maximum) ) or ( $maximum < 0 ) ) {
    die "Maximum ortholog number must be a positive integer, not: $maximum\n";
}

# To be used as filehandle for the input(s):
my $INFILE;

foreach my $infile (@infiles) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INFILE>) { 
        chomp $input;

        # Sample input line:  
        # ORTHOMCL19007(2 genes,2 taxa: WBGene00031156(briggsae) WBGene00005532(elegans)

        if ( ( $input =~ /\S/xms ) 
             and ( $input !~ / \A 
                               [A-Za-z]+\d+ 
                               \( \d+ \s+ genes,\d+ \s+ taxa\) : 
                               (?: \s+ [^\s\(\)]+ \( [^\s\(\)]+ \) )+
                               \s* \z /xms ) ) { 
            die "Obviously not an OrthoMCL or OrthoFinder text line!\n";
        }
        else { 
            while ( $input =~ / (\S+) \( $taxon \) /gxms ) { 
                my $gene = $1;
                my $orthomcl = $input;
                if ($summary) { 
                    $orthomcl = summ_omcl($gene,$orthomcl);
                }
                # Native OrthoMCL file format has a '\t' on the line, which is pretty, but which kills *.tsv registry; so s/ it out.
                if (! $summary ) { 
                    $orthomcl =~ s/\t/ /g;
                }
                # Truncate reallllllly long lines, if $maximum is set.
                if ($maximum) { 
                     my $omcl_length = length $orthomcl;
                     my $excess_omcl = $omcl_length - $maximum;

                     $orthomcl =  substr("$orthomcl", 0, $maximum);

                     if ( $excess_omcl >= 1 ) { 
                         my $excess_text = " [etc. -- $excess_omcl characters omitted]";
                         $orthomcl .= $excess_text;
                     }
                }
                $data_ref->{'gene'}->{$gene}->{'omcl'}->{$orthomcl} = 1;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile. $!\n";
}

my $header = "Gene\t";
if ($summary) { 
    $header .= "OFind_Summary\n";
}
if (! $summary) {
    $header .= "OFind_Full\n";
}

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @omcls = ();
    if ( exists $data_ref->{'gene'}->{$gene1}->{'omcl'} ) { 
        @omcls = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'omcl'} };
    }
    my $omcl_text = join '; ', @omcls;
    print $header if $header;
    $header = q{};
    print "$gene1\t$omcl_text\n";
}

sub summ_omcl { 
    my $_gene     = $_[0];
    my $_orthomcl = $_[1];
    if ( $_orthomcl !~ / \A [A-Za-z]+\d+ \( \d+ \s+ genes,\d+ \s+ taxa\) : (?: \s+ [^\s\(\)]+ \( [^\s\(\)]+ \) )+ \s* \z /xms ) {
        die "Subroutine summ_omcl can't parse input: $_orthomcl\n";
    }
    my $_omcl_name = q{};
    if ( $_orthomcl =~ / \A ( [A-Za-z]+\d+ \( \d+ \s+ genes,\d+ \s+ taxa\) ) : ( (?: \s+ [^\s\(\)]+ \( [^\s\(\)]+ \) )+ ) \s* \z /xms ) {
        $_omcl_name       = $1;
        my $_omcl_members = $2;

        while ( $_omcl_members =~ / \S+ \( ( \S+ ) \) /gxms ) {
            my $_taxon = $1;
            $data_ref->{'gene'}->{$_gene}->{'_omcl_name'}->{$_omcl_name}->{'_taxon'}->{$_taxon} += 1;
        }
    }
    my @_taxa = sort keys %{ $data_ref->{'gene'}->{$_gene}->{'_omcl_name'}->{$_omcl_name}->{'_taxon'} };
    my @_omcl_taxon_counts = ();
    foreach my $_taxon1 (@_taxa) { 
        my $_taxon_no    = $data_ref->{'gene'}->{$_gene}->{'_omcl_name'}->{$_omcl_name}->{'_taxon'}->{$_taxon1};
        my $_taxon_count = "$_taxon1 ($_taxon_no g.)";
        push @_omcl_taxon_counts, $_taxon_count;
    }
    my $_omcl_taxon_text = join '; ', @_omcl_taxon_counts;
    my $_omcl_summary = $_omcl_name . ': ' . $_omcl_taxon_text;
    return $_omcl_summary;
}

