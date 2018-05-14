#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);
use Scalar::Util qw(looks_like_number);

my $data_ref;

my $L3_data    = q{};
my $L4_data    = q{};
my $N2_data    = q{};
my $nhr67_data = q{};

my %important  = ();

my $help;

GetOptions ( 'L3=s'  => \$L3_data,
             'L4=s'  => \$L4_data,
             'N2=s'  => \$N2_data,
             'nhr=s' => \$nhr67_data,
             'help'  => \$help,      );

my %infiles     = ( L3  => "$L3_data",
                    L4  => "$L4_data",
                    N2  => "$N2_data",
                    nhr => "$nhr67_data", );

my $filecount = grep { $_ =~ /\S/ } values %infiles;

if ($help or ( $filecount != 4 ) ) { 
    die "Format: single_LC_table_21jun2010.pl\n",
        "    --L3   [wild-type L3 RPKM data file]\n",
        "    --L4   [wild-type L4 RPKM data file]\n",
        "    --N2   [N2 larval RPKM data file]\n",
        "    --nhr  [nhr-67(RNAi) L4 RPKM data file]\n",
        ;
}

foreach my $expr_type (sort keys %infiles) { 
    my $expr_data = $infiles{$expr_type};
    open my $EXPR_DATA, '<', $expr_data or die "Can't open data file $expr_data for expression type $expr_type: $!";
    my $gene = q{};
    my $rpkm = q{};
    while (my $input = <$EXPR_DATA>) { 
        chomp $input;
        if ( $input =~ / \A (\S+) \t [^\t]* \t (\d+\.\d+) /xms) { 
            $gene = $1;
            $rpkm = $2;
            if (! ( looks_like_number $rpkm ) ) {
                die "RPKM $rpkm does not look numerical.\n";
            }
            # Avoid having real data overwritten by bogus '0.00' data, from split CDSes!
            if ( exists $data_ref->{$gene}->{$expr_type}->{'RPKM'} ) { 
                my $prev_rpkm = $data_ref->{$gene}->{$expr_type}->{'RPKM'};
                if ( $rpkm > $prev_rpkm ) { 
                    $data_ref->{$gene}->{$expr_type}->{'RPKM'} = $rpkm;
                }
            }
            if (! exists $data_ref->{$gene}->{$expr_type}->{'RPKM'} ) {
                $data_ref->{$gene}->{$expr_type}->{'RPKM'} = $rpkm;
            }
            # This script is being used specifically to get only wild-type LC active genes.
            if ( ( ( $expr_type eq 'L3' ) or ( $expr_type eq 'L4') ) and ( $rpkm > 0 ) ) { 
                $important{$gene} = 1;
            }
        }
        else { 
            if ( $input ne "#gene\tlen_kb\tRPKM\tmulti/all" ) {
                warn "From data file $expr_data for expression",
                     " type $expr_type, can't parse input line:\n",
                     "$input\n",
                  ;
            }
        }
    }
    close $EXPR_DATA 
        or die "Can't close filehandle for data file $expr_data for expression type $expr_type: $!";
}

foreach my $imp_gene (sort keys %important) { 
    # Enforce *some* value; specifically, '0.00' for anything not recorded.
    foreach my $expr_type2 (sort keys %infiles) { 
        if (! exists $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM'} ) { 
            $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM'} = 0.00;
        }
    }

    # Compute the maximum wild-type to N2 expression ratio:
    my $max_rpkm = max( $data_ref->{$imp_gene}->{'L3'}->{'RPKM'}, 
                        $data_ref->{$imp_gene}->{'L4'}->{'RPKM'}  );
    my $bg_rpkm = $data_ref->{$imp_gene}->{'N2'}->{'RPKM'};
    $data_ref->{$imp_gene}->{'max_RPKM'} = $max_rpkm;

    # Give the background N2 RPKM a nominal minimum of 0.01 to avoid division by zero.
    if ( $bg_rpkm == 0.00 ) { 
        $bg_rpkm = 0.01;
    }

    # Finally, get the ratio, in human-readable form:
    my $ratio = $max_rpkm / $bg_rpkm;
    $ratio = sprintf "%.2f", $ratio;
    my $rpkm_fields =  $data_ref->{$imp_gene}->{'L3'}->{'RPKM'}
                      . "\t"
                      . $data_ref->{$imp_gene}->{'L4'}->{'RPKM'}
                      . "\t"
                      . $data_ref->{$imp_gene}->{'nhr'}->{'RPKM'}
                      . "\t" 
                      . $data_ref->{$imp_gene}->{'N2'}->{'RPKM'}
                      ;
    $data_ref->{'final'}->{$imp_gene}->{'ratio'} = $ratio;
    $data_ref->{'final'}->{$imp_gene}->{'RPKMs'} = $rpkm_fields;
}

my @final_genes = sort { $data_ref->{'final'}->{$b}->{'ratio'} 
                         <=> $data_ref->{'final'}->{$a}->{'ratio'} } 
                  sort { $data_ref->{$b}->{'max_RPKM'} 
                         <=> $data_ref->{$a}->{'max_RPKM'} } 
                  keys %{ $data_ref->{'final'} };

foreach my $final_gene (@final_genes) { 
    print $final_gene,
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'ratio'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'RPKMs'},
          "\n",
          ;
}

