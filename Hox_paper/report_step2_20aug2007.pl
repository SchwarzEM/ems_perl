#!/usr/bin/perl

# report_step2_v02.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/19/2007
# Purpose: generate second stage of Tables 1, 2 ("report") on ~8/19/2007.

use strict;
use warnings;

unless ($#ARGV == 2) { 
    die "Format: ./report_step2.pl ",
        "[gene names] [orthomcl] ",
        "[report step 1]\n"
        ;
}

# Muster relevant filenames.
my $gene_names   = $ARGV[0];
my $orthomcl     = $ARGV[1];
my $report_step1 = $ARGV[2];

# Map WBGene\d+ to interesting C. elegans names.
my %gene2name    = ();

# Store orthology data for later reworking.
my %gene2ortho   = ();

# And link prot-IDs to their full-formatted ortho. descs.
my %gene2o_desc  = ();

open my $GENENAMES, $gene_names 
    or die "Can't open gene name file $gene_names: $!";
while (my $input1 = <$GENENAMES>) { 
    chomp $input1;
    if ($input1 =~ /\A (WBGene\d+)     # $1
                      .* \s 
                      (\S+-\S+)        # $2: '-' == CGC name 
                      \s* /xms) {
        $gene2name{$1} = $2;
    }
    elsif ($input1 =~ /\A (WBGene\d+)  # $1 
                      .* \s 
                      (\S.+\S)         # $2: '.' == CDS name
                      \s* /xms) {
        $gene2name{$1} = $2;
    }
}
close $GENENAMES;

open my $ORTHOMCL, $orthomcl 
    or die "Can't open OrthoMCL report $orthomcl: $!";
while (my $input2 = <$ORTHOMCL>) { 
    chomp $input2;

    # If this line from a vast OrthoMCL output actually has CB5161/PS1010...
    if ($input2 =~ / cb5161 | ps1010 /xms) { 

        # Rename WBGene\d+ IDs to best human-readable names.
        if ($input2 =~ /WBGene/xms) { 
            $input2 =~ s/(WBGene\d+)/$gene2name{$1}/g;
        }

        # Clip off header information.
        $input2 =~ s/ \A ORTHOMCL\d+ 
                      \( 
                      \d+\sgenes,\d+\staxa
                      \):
                      \s+
                    //xms;

        # Now, extract most salient details of orthologs -- or try to!
        my @o_genes = split /\s/, $input2;
        foreach my $o_gene (@o_genes) { 
            if ( $o_gene =~ /\A 
                             (\S+) 
                             \(
                             (cb5161|ps1010)
                             \)
                             /xms ) { 
                @{ $gene2ortho{$1} } = @o_genes;
            }
        }
    }
}
close $ORTHOMCL;

foreach my $o_gene (sort keys %gene2ortho) { 

    # Initialize: counts of genes/spec.; key orthologs; final desc. text.
    my %sp_count          = ();
    my %interesting_orths = ();
    my $full_output       = q{};

    foreach my $ortholog (sort @{ $gene2ortho{$o_gene} }) { 
        if ($o_gene ne $gene2ortho{$o_gene}) { 
            my $species = q{};
            if ( $ortholog =~ / 
                                ( elegans  | 
                                  briggsae | 
                                  remanei  | 
                                  cb5161   | 
                                  ps1010 )
                               /xms ) { 
                $species = $1;
                $sp_count{$species}++;
            }
            if ( $ortholog =~ /\A 
                               (\S+)
                               \(
                               (elegans)   # Usually only 'elegans'...
                               \)
                               /xms ) { 
                push @{ $interesting_orths{"$2"} }, $1;
            }
        }
    }

    # List interesting orthologs:
    foreach my $species (sort keys %interesting_orths) {

        # "Interesting" means elegans, of course:
        if ( $interesting_orths{$species} ) {
            my $cool_output = q{};
            if (@{ $interesting_orths{$species} } >= 2) { 
                $interesting_orths{$species}->[-1] 
                    = "and " . $interesting_orths{$species}->[-1];
            }
            if ( @{ $interesting_orths{$species} } >= 3 ) {
                $cool_output = join ", ", @{ $interesting_orths{$species} };
            }
            if ( @{ $interesting_orths{$species} } == 2 ) {
                $cool_output = join " ", @{ $interesting_orths{$species} };
            }
            if ( @{ $interesting_orths{$species} } == 1 ) {
                $cool_output = $interesting_orths{$species}->[0] ;
            }

            # Start building up the ortholog description text.
            $full_output .= "$cool_output";
        }

        # Note the elegans orthologs that are highly unambiguous:
        if (  ( $sp_count{"elegans"} and ( $sp_count{"elegans"} == 1 ) ) 
              and 
              (   ( $sp_count{"briggsae"} and ( $sp_count{"briggsae"} == 1 ) )
                  or 
                  ( $sp_count{"remanei"}  and ( $sp_count{"remanei"}  == 1 ) )
              )
           ) { 

            # Add more ortholog text:
            $full_output .= ' [*]';
        }
        # Add more ortholog text:
        $full_output .= " ";
    }

    # List general orthologs:
    my @gen_list = ();
    foreach my $species qw(elegans briggsae remanei cb5161 ps1010) { 
        if ($sp_count{$species}) { 
            push @gen_list, "$sp_count{$species} $species";
        }
    }
    my $gen_output = join ", ", @gen_list;

    # Add more ortholog text.
    $full_output .= "($gen_output).";  # No \n at end.

    # Link each gene to its ortholog description.
    $gene2o_desc{$o_gene} = $full_output;
}

open my $REPORT1, $report_step1
    or die "Can't open preliminary report $report_step1: $!";
while (my $input = <$REPORT1>) { 
    chomp $input;
    print $input;
    if ($input =~ /\A (\S+) /xms) { 
        my $gene = $1;
        if ( $gene2o_desc{$gene} ) { 
            print "\t$gene2o_desc{$gene}";
        }
    }
    print "\n";
}
# Give blank line at bottom.
print "\n";
close $REPORT1;

