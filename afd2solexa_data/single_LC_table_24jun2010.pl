#!/usr/bin/env perl

# single_LC_table_24jun2010.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/24/2010.
# Purpose: make a TSV file of four RPKM values and their ratios (absolute and log[10]).

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

    # Get concise strings for individual RPKM data:
    my $L3_rpkm    = $data_ref->{$imp_gene}->{'L3'}->{'RPKM'};
    my $L4_rpkm    = $data_ref->{$imp_gene}->{'L4'}->{'RPKM'};
    my $nhr67_rpkm = $data_ref->{$imp_gene}->{'nhr'}->{'RPKM'};
    my $N2_rpkm    = $data_ref->{$imp_gene}->{'N2'}->{'RPKM'};

    # Get a concise string variable for maximum wild-type LC RPKM:
    my $maxLC_rpkm = max( $data_ref->{$imp_gene}->{'L3'}->{'RPKM'}, 
                        $data_ref->{$imp_gene}->{'L4'}->{'RPKM'}  );
    # Use this as part of sorting the output, later:
    $data_ref->{$imp_gene}->{'maxLC_RPKM'} = $maxLC_rpkm;

    # Archive a single printable text-line segment of RPKM values:
    my $basic_rpkms = $L3_rpkm
                    . "\t"
                    . $L4_rpkm
                    . "\t"
                    . $nhr67_rpkm
                    . "\t"
                    . $N2_rpkm
                    . "\t"
                    . $maxLC_rpkm
                    ;

    $data_ref->{'final'}->{$imp_gene}->{'RPKMs'} = $basic_rpkms;

    # Get log[10] versions of these worked up and stored as text-lineseg.
    # Can use safely_divide_RPKMs, since getting the log of X is equivalent to getting the log of X/1!
    #     In other words, the job I might want to do with safely_take_log is a subset of what safely_divide_RPKMs does.
    #     So, just invoke that with a simple hack...

    my $L3_log10    = safely_divide_RPKMs($L3_rpkm,    1, 0.01, 'log', 10);
    $L3_log10       = sprintf "%.3f", $L3_log10;

    my $L4_log10    = safely_divide_RPKMs($L4_rpkm,    1, 0.01, 'log', 10);
    $L4_log10       = sprintf "%.3f", $L4_log10;

    my $nhr_log10   = safely_divide_RPKMs($nhr67_rpkm, 1, 0.01, 'log', 10);
    $nhr_log10      = sprintf "%.3f", $nhr_log10;

    my $N2_log10    = safely_divide_RPKMs($N2_rpkm,    1, 0.01, 'log', 10);
    $N2_log10       = sprintf "%.3f", $N2_log10;

    my $maxLC_log10 = safely_divide_RPKMs($maxLC_rpkm, 1, 0.01, 'log', 10);
    $maxLC_log10    = sprintf "%.3f", $maxLC_log10;

    my $basic_rpkm_log10s = $L3_log10
                          . "\t"
                          . $L4_log10
                          . "\t"
                          . $nhr_log10
                          . "\t"
                          . $N2_log10
                          . "\t"
                          . $maxLC_log10
                          ;

    $data_ref->{'final'}->{$imp_gene}->{'RPKM_log10s'} = $basic_rpkm_log10s;

    # Get the max./N2 ratio data, in human-readable plain and log[10] form:
    my $maxLC_N2_ratio = safely_divide_RPKMs($maxLC_rpkm, $N2_rpkm, 0.01);
    $maxLC_N2_ratio = sprintf "%.2f", $maxLC_N2_ratio;

    my $maxLC_N2_log10_ratio = safely_divide_RPKMs($maxLC_rpkm, $N2_rpkm, 0.01, 'log', 10);
    $maxLC_N2_log10_ratio = sprintf "%.3f", $maxLC_N2_log10_ratio;

    # Store these data individually, so that I can later use them to sort the output:
    $data_ref->{'final'}->{$imp_gene}->{'maxLC_vs_N2_ratio'} = $maxLC_N2_ratio;
    $data_ref->{'final'}->{$imp_gene}->{'maxLC_vs_N2_log10_ratio'} = $maxLC_N2_log10_ratio;

    # Store this as a single printable text-line segment:
    $data_ref->{'final'}->{$imp_gene}->{'maxLC_vs_N2_data'} 
        = $maxLC_N2_ratio
          . "\t"
          . $maxLC_N2_log10_ratio
          ;

    # Get the wt L4/wt L3 ratio data, and store as printable text-line segment:
    my $L4_L3_ratio = safely_divide_RPKMs($L4_rpkm, $L3_rpkm, 0.01);
    $L4_L3_ratio = sprintf "%.2f", $L4_L3_ratio;
    my $L4_L3_log10_ratio = safely_divide_RPKMs($L4_rpkm, $L3_rpkm, 0.01, 'log', 10);
    $L4_L3_log10_ratio = sprintf "%.3f", $L4_L3_log10_ratio;
    $data_ref->{'final'}->{$imp_gene}->{'L4_vs_L3_data'} 
        = $L4_L3_ratio
          . "\t"
          . $L4_L3_log10_ratio
          ;

    # Get the nhr-67 L4/wt L3 ratio data, and store as printable text-line segment:
    my $nhr67_L3_ratio = safely_divide_RPKMs($nhr67_rpkm, $L3_rpkm, 0.01);
    $nhr67_L3_ratio = sprintf "%.2f", $nhr67_L3_ratio;
    my $nhr67_L3_log10_ratio = safely_divide_RPKMs($nhr67_rpkm, $L3_rpkm, 0.01, 'log', 10);
    $nhr67_L3_log10_ratio = sprintf "%.3f", $nhr67_L3_log10_ratio;
    $data_ref->{'final'}->{$imp_gene}->{'nhr67_vs_L3_data'}
        = $nhr67_L3_ratio
          . "\t"
          . $nhr67_L3_log10_ratio
          ;

    # Get the nhr-67 L4/wt L4 ratio data, and store as printable text-line segment:
    my $nhr67_L4_ratio = safely_divide_RPKMs($nhr67_rpkm, $L4_rpkm, 0.01);
    $nhr67_L4_ratio = sprintf "%.2f", $nhr67_L4_ratio;
    my $nhr67_L4_log10_ratio = safely_divide_RPKMs($nhr67_rpkm, $L4_rpkm, 0.01, 'log', 10);
    $nhr67_L4_log10_ratio = sprintf "%.3f", $nhr67_L4_log10_ratio;
    $data_ref->{'final'}->{$imp_gene}->{'nhr67_vs_L4_data'}
        = $nhr67_L4_ratio
          . "\t"
          . $nhr67_L4_log10_ratio
          ;           



}

my @final_genes = sort { $data_ref->{'final'}->{$b}->{'maxLC_vs_N2_ratio'} 
                         <=> $data_ref->{'final'}->{$a}->{'maxLC_vs_N2_ratio'} } 
                  sort { $data_ref->{'final'}->{$b}->{'maxLC_vs_N2_log10_ratio'}
                         <=> $data_ref->{'final'}->{$a}->{'maxLC_vs_N2_log10_ratio'} }
                  sort { $data_ref->{$b}->{'maxLC_RPKM'} 
                         <=> $data_ref->{$a}->{'maxLC_RPKM'} } 
                  keys %{ $data_ref->{'final'} };

foreach my $final_gene (@final_genes) { 
    print $final_gene,
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'RPKMs'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'RPKM_log10s'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'maxLC_vs_N2_data'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'L4_vs_L3_data'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'nhr67_vs_L3_data'},
          "\t",
          $data_ref->{'final'}->{$final_gene}->{'nhr67_vs_L4_data'},
          "\n",
          ;
}

sub safely_divide_RPKMs { 
    my $_numerator   = $_[0];
    my $_denominator = $_[1];
    my $_safezero    = $_[2];
    my $_take_log    = $_[3];   # Any string (e.g., 'log') is a Boolean 'yes'; default is 'no'.
    my $_log_base    = $_[4];

    # Safeguard the safezero: 
    if ( $_safezero == 0 ) { 
        die "When dividing $_numerator by $_denominator,",
            " cannot safely use $_safezero as a substitute for zero.\n",
            ;
    }

    # Refuse to use zero or negative numbers as a logarithmic base:
    if ( $_take_log and ( $_log_base <= 0 ) ) { 
        die "When getting log of $_numerator/$_denominator,",
        " cannot safely use $_log_base as a log base.\n",
        ;
    }

    # No division by zero:
    if ( $_denominator == 0 ) { 
        $_denominator = $_safezero;
    }

    # No logarithms of zero:
    if ($_take_log and ( $_numerator == 0 ) ) {
        $_numerator = $_safezero;
    }

    # Default without a log base: compute a simple ratio.
    my $_ratio = ($_numerator / $_denominator);

    # Given a positive log base: get its log (usually will be log[10]).
    if ($_take_log) { 
        $_ratio = ( log($_ratio) / log($_log_base) );
    }

    # Return the simple ratio or its log.
    return $_ratio;
}

