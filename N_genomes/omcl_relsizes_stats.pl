#!/usr/bin/env perl

# omcl_relsizes_stats.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/19/2009.
# Purpose: for one or more OrthoMCLs, overall statistics of ratios of each sp. longest protein to elegans's longest protein.

use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;

my %proteomes      = ();
my %prot2len       = ();
my $orthomcl_input = q{};

# Standard order of recording/reading species data:
my @spp = qw( elegans briggsae remanei brenneri japonica );

# Reference to species-specific array of length values:
my $spp2len_ref;

GetOptions (  'ele:s' => \$proteomes{'elegans'},
              'bri:s' => \$proteomes{'briggsae'},
              'rem:s' => \$proteomes{'remanei'},
              'bre:s' => \$proteomes{'brenneri'},
              'jap:s' => \$proteomes{'japonica'},

              'omcl:s' => \$orthomcl_input,         
           );

### Parse input files: ###

my $proteome_count = keys %proteomes;

unless (     ( $proteome_count >= 2  ) 
         and ( $proteomes{'elegans'} )
         and ( $orthomcl_input       ) ) { 
    die_loudly();
}

### Read and store protein lengths for each proteome. ### 

foreach my $species (@spp) { 
    if ($proteomes{$species}) { 
        open my $PROTEOME, '<', $proteomes{$species} 
            or die "Cannot open proteome $proteomes{$species}: $!";
        my $cds = q{};
        while (my $input = <$PROTEOME>) {
            chomp $input;
            if ($input =~ /\A > (\S+) \s* /xms) {
                $cds    = $1;
                # Ensure OrthoMCL suffix on $cds name:
                if ( $cds !~ / \( $species \) \z  /xms ) { 
                    $cds .= "($species)";
                }
                # And then start the residue count.
                $prot2len{$cds} = 0;
            }
            elsif ( $input =~ /\S/xms ) {
                # Absolutely minimal filtering:
                $input =~ s/\s//g;
                # Increment residue count:
                $prot2len{$cds} += length($input);
            }
        }
        close $PROTEOME
            or die "Can't close filehandle to proteome",
                   " $proteomes{$species}: $!",
                   ;
    }
}

### Get max-res counts and ratios vs. elegans: ###

# Sample input:
# ORTHOMCL5940(5 genes,5 taxa):	 CBG28141(briggsae) CBN08237(brenneri) CJA27759(japonica) CRE00189(remanei) T13C5.8(elegans)

open my $ORTHO_INPUT, '<', $orthomcl_input 
    or die "Cannot open N-species OrthoMCL output $orthomcl_input: $!";

while (my $input = <$ORTHO_INPUT>) { 
    chomp $input;
    if ($input =~ / \A
                    (ORTHOMCL\d+)                 # $1 -> $ortho_grp
                    \(\d+\sgenes,\d+\staxa\)      # just (punctuation)
                    : \s+ 
                    (.+)                          # $2 -> $oprots_line
                    \s* 
                    \z 
                  /xms) { 
        my $ortho_grp = $1;
        my $oprots_line = $2;
        my %species2max = ();
        foreach my $sp1 (@spp) { 
            $species2max{$sp1} = 0;
        }
        my @orthoprots = split /\s+/, $oprots_line;
        foreach my $o_prot (@orthoprots) { 
            if ( $o_prot !~ / .+ \( \w+ \) /xms ) {
                die "Can't parse OrthoMCL protein $o_prot\n";
            }
            if ( $o_prot =~ / .+ \( (\w+) \) /xms ) { 
                my $species_tag = $1;
                if (! $prot2len{$o_prot} ) { 
                    die "Failed to record length of protein $o_prot!\n";
                }
                my $prot_size = $prot2len{$o_prot};
                if ( $prot_size > $species2max{$species_tag} ) { 
                    $species2max{$species_tag} = $prot_size;
                }
            }
        }
        my @length_data = ();
        push @length_data, $ortho_grp;
        if ( $species2max{'elegans'} >= 1 ) { 
            foreach my $sp2 (@spp) { 
                my $ratio = ($species2max{$sp2} / $species2max{'elegans'} );
                $ratio = sprintf "%.3f", $ratio;
                push @{ $spp2len_ref->{$sp2}->{'ratio'} }, $ratio;
                push @{ $spp2len_ref->{$sp2}->{'aa_count'} }, 
                    $species2max{$sp2};
            }
        }
    }
}
close $ORTHO_INPUT 
    or die "Cannot close filehandle to $orthomcl_input: $!";

foreach my $sp3 (@spp) { 
    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data( @{ $spp2len_ref->{$sp3}->{'ratio'} } );
    my $stat2 = Statistics::Descriptive::Full->new();
    $stat2->add_data( @{ $spp2len_ref->{$sp3}->{'aa_count'} } );

    my $omcls = $stat1->count();  # Total number of contigs.
    $omcls = commify($omcls);

    my $total_aa = $stat2->sum();
    $total_aa = commify($total_aa);

    my $mean    = $stat1->mean();
    $mean       = sprintf("%.3f", $mean);

    my $std_dev = $stat1->standard_deviation();
    $std_dev    = sprintf("%.3f", $std_dev);  

    my $min     = $stat1->min();
    my $max     = $stat1->max();
    my $median  = $stat1->median();

    print "\n";
    print "Species $sp3: $omcls OrthoMCL groups with $total_aa total a.a.\n";
    print "Mean:   $mean; std. dev. $std_dev; min. $min; max. $max\n";
    print "Median: $median\n";
    print "\n";
}

### Subroutines: ###

sub die_loudly {
    die 'Format: ./edit_3worm_orthomcl.pl                 \ ', "\n",
        '    # Need elegans and at least other proteome:  \ ', "\n",
        '    --ele  [elegans proteome, REQUIRED]          \ ', "\n",
        '    --bri  [briggsae proteome]                   \ ', "\n",
        '    --rem  [remanei proteome]                    \ ', "\n",
        '    --bre  [brenneri proteome]                   \ ', "\n",
        '    --jap  [japonica proteome]                   \ ', "\n",
        '    # Need a single OrthoMCL output:             \ ', "\n",
        '    --ocml [orthomcl.out]                        \ ', "\n",
        ;
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

