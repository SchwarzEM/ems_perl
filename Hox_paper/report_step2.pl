#!/usr/bin/perl

# report_step2.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/10/2007
# Purpose: generate second stage of Tables 1, 2 ("report") on 12/10/2007.
# Revised from earlier version on 8/xx/2007.

use strict;
use warnings;

unless ($#ARGV == 2) { 
    die "Format: ./report_step2.pl ",
        "[gene names] [orthomcl] ",
        "[report step 1]\n"
        ;
}

# Table for replacing old with new CB5161/PS1010 clone/gene names:
my %renamed = ( 
    'CB5161_anon1.tfa'         => "Cbre_JD01",
    'CB5161_anon2.tfa'         => "Cbre_JD02",
    'CB5161_anon3A.tfa'        => "Cbre_JD03",
    'CB5161_anon3B.tfa'        => "Cbre_JD04",
    'CB5161_anon4A.tfa'        => "Cbre_JD05",
    'CB5161_anon4B.tfa'        => "Cbre_JD06",
    'CB5161_ceh-13.lin-39.tfa' => "Cbre_JD07",
    'CB5161_ceh-8.tfa'         => "Cbre_JD08",
    'CB5161_col-14.tfa'        => "Cbre_JD09",
    'CB5161_egl-44.tfa'        => "Cbre_JD10",
    'CB5161_egl-46.tfa'        => "Cbre_JD11",
    'CB5161_gpa-6.tfa'         => "Cbre_JD12",
    'CB5161_lin-11.tfa'        => "Cbre_JD13",
    'CB5161_lin-3.tfa'         => "Cbre_JD14",
    'CB5161_ndescA.tfa'        => "Cbre_JD15",
    'CB5161_ndescB.tfa'        => "Cbre_JD16",
    'CB5161_nlp-8.tfa'         => "Cbre_JD17",
    'CB5161_rde-1A.tfa'        => "Cbre_JD18",
    'CB5161_rde-1B.tfa'        => "Cbre_JD19",
    'CB5161_ref-1.tfa'         => "Cbre_JD20",
    'CB5161_srw-2.tfa'         => "Cbre_JD21",
    'CB5161_weird.tfa'         => "Cbre_JD22",

    'PS1010_anonymous.tfa'     => "Csp3_JD01",
    'PS1010_ceh-13.lin-39.tfa' => "Csp3_JD02",
    'PS1010_egl-30.tfa'        => "Csp3_JD03",
    'PS1010_egl-5.mab-5.tfa'   => "Csp3_JD04",
    'PS1010_lin-11.tfa'        => "Csp3_JD05",
    'PS1010_lin-3.tfa'         => "Csp3_JD06",
    'PS1010_php-3.tfa'         => "Csp3_JD07",
);

# Invoke this as needed to update our names, best right before print:
sub rename_our_genes { 
    my $_input = $_[0];
    my @old_names = keys %renamed;
    foreach my $old_name (@old_names) { 
        $_input =~ s/$old_name/$renamed{$old_name}/g;
    }
    return $_input;
}

sub rename_species { 
    my $_input = $_[0];
    $_input =~ s/cb5161/brenneri/g;
    # Could rename 'ps1010' if needed.
    return $_input;
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
            $input2 =~ s/(WBGene\d+)/$1|$gene2name{$1}/g;
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
            if ( @{ $interesting_orths{$species} } >= 2 ) { 
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
    my $output = $input;
    if ($input =~ /\A (\S+) /xms) { 
        my $gene = $1;
        if ( $gene2o_desc{$gene} ) { 
            $output .= "\t$gene2o_desc{$gene}";
        }
    }
    $output = rename_our_genes($output);
    $output = rename_species($output);
    print "$output\n";
}
# Give blank line at bottom.
print "\n";
close $REPORT1;

