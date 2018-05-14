#!/usr/bin/env perl

# single_LC_table_12oct2010.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/12/2010.
# Purpose: make *.tsv of linker RPKM values and ratios (abs. and log[10]); simple and non-multi. (uniq+splice) RPKMs; optional raw, unrounded values; optional header text.

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

my %important_gene      = ();
my @important_data      = ();
my %important_expr_type = ();
my %ok_expr_type = ( L3  => 1,
                     L4  => 1,
                     N2  => 1,
                     nhr => 1, );

my $help;
my $no_rounding;
my $header_text;
my $auxillary;

GetOptions ( 'L3=s'     => \$L3_data,
             'L4=s'     => \$L4_data,
             'N2=s'     => \$N2_data,
             'nhr=s'    => \$nhr67_data,
             'imp:s{,}' => \@important_data,
             'header'   => \$header_text,
             'raw'      => \$no_rounding,
             'aux'      => \$auxillary,
             'help'     => \$help,      );

# This default gets re-used later, when deciding what maximum values are interesting to have in ratios with N2
#     (as opposed to which maximum values are worth using in the sorting of the genes).
if (! @important_data) { 
    @important_data = qw(L3 L4);
}
foreach my $imp_dat (sort @important_data) { 
    if (! exists $ok_expr_type{$imp_dat} ) { 
        die "The '--imp' argument only takes one or more of: 'L3', 'L4', 'N2', or 'nhr'.\n";
    }
    $important_expr_type{$imp_dat} = 1;
}

my %infiles     = ( L3  => "$L3_data",
                    L4  => "$L4_data",
                    N2  => "$N2_data",
                    nhr => "$nhr67_data", );

my $filecount = grep { $_ =~ /\S/ } values %infiles;

if ($help or ( $filecount != 4 ) ) { 
    die "\n",
        "Format: single_LC_table_12oct2010.pl\n",
        "    --L3      [wild-type L3 RPKM data file]\n",
        "    --L4      [wild-type L4 RPKM data file]\n",
        "    --N2      [N2 larval RPKM data file]\n",
        "    --nhr     [nhr-67(RNAi) L4 RPKM data file]\n",
        "    --imp     [important expression types; defaults to 'L3' and 'L4'; can pick any comb. of 'L3', 'L4', 'N2', or 'nhr']\n",
        "    --header  [optionally, insert header line at top of TSV file to explain columns]\n",
        "    --raw     [do not round to 0.01; useful for fine-resolution data plotting]\n",
        "    --aux     [include auxillary data about highest RPKMs of important expression types]\n",
        "    --help    [display these options]\n",
        "\n",
        "    Purpose:   make a TSV file of four RPKM values and their ratios (absolute and log[10]),\n",
        "               with both simple RPKMs and non-multiread (uniq + splice) RPKMs;\n",
        "               optionally export raw, unrounded values rather than values rounded to 0.01.\n",
        "\n",
        ;
}

foreach my $expr_type (sort keys %infiles) { 
    my $expr_data = $infiles{$expr_type};
    open my $EXPR_DATA, '<', $expr_data or die "Can't open data file $expr_data for expression type $expr_type: $!";
    my $gene           = q{};
    my $rpkm           = q{};
    my $mult_frac      = q{};
    my $rpkm_non_multi = q{};   # It includes splices, which in principle *could* map to >= 2 sites, so it's not "unique-only".
    while (my $input = <$EXPR_DATA>) { 
        chomp $input;

# Sample input:
#   #gene	len_kb	RPKM	multi/all
#   WBGene00009119|F25H2.5	0.818	98437.54	0.00

        if ( $input =~ / \A (\S+) \t [^\t]* \t (\d+\.\d+) \t (\d.\d{2}) /xms) { 
            $gene      = $1;
            $rpkm      = $2;
            $mult_frac = $3;
            if ( (! ( looks_like_number $rpkm ) ) or (! ( looks_like_number $mult_frac ) ) ) {
                die "In $input: RPKM $rpkm or multi/all ratio $mult_frac (or both) does not look numerical.\n";
            }
            if ( ( $mult_frac < 0 ) or ( $mult_frac > 1 ) ) { 
                die "In $input: insane multi/all ratio $mult_frac.\n";
            }
            $rpkm_non_multi = ($rpkm * ( 1.00 - $mult_frac ));

            # Avoid having real data overwritten by bogus '0.00' data, from split CDSes!
            if ( exists $data_ref->{$gene}->{$expr_type}->{'RPKM'} ) { 
                my $prev_rpkm = $data_ref->{$gene}->{$expr_type}->{'RPKM'};
                # Also, keep raising the $rpkm until maximized; pass over submaximal values.
                if ( $rpkm > $prev_rpkm ) { 
                    $data_ref->{$gene}->{$expr_type}->{'RPKM'} = $rpkm;
                }
            }
            if (! exists $data_ref->{$gene}->{$expr_type}->{'RPKM'} ) {
                $data_ref->{$gene}->{$expr_type}->{'RPKM'} = $rpkm;
            }

            # Do likewise with non-multiread RPKMs:
            if ( exists $data_ref->{$gene}->{$expr_type}->{'RPKM_non_multi'} ) {
                my $prev_rpkm_non_multi = $data_ref->{$gene}->{$expr_type}->{'RPKM_non_multi'};
                # Excelsior!
                if ( $rpkm_non_multi > $prev_rpkm_non_multi ) {
                    $data_ref->{$gene}->{$expr_type}->{'RPKM_non_multi'} = $rpkm_non_multi;
                }
            }
            if (! exists $data_ref->{$gene}->{$expr_type}->{'RPKM_non_multi'} ) {
                $data_ref->{$gene}->{$expr_type}->{'RPKM_non_multi'} = $rpkm_non_multi;;
            }

            # By default, this script is being used specifically to get only wild-type LC active genes ('L3' or 'L4').
            # However, the user can specify other sets.
            if ( ( exists $important_expr_type{$expr_type} ) and ( $rpkm > 0 ) ) {
                $important_gene{$gene} = 1;
            }
        }
        else { 
            if ( $input ne "#gene\tlen_kb\tRPKM\tmulti\/all" ) {
                # Don't die, presumably because there may be sporadic bad lines in mostly good data.
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

foreach my $imp_gene (sort keys %important_gene) { 
    # Enforce *some* value; specifically, '0.00' for anything not recorded.
    foreach my $expr_type2 (sort keys %infiles) { 
        if (! exists $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM'} ) { 
            $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM'} = 0.00;
        }
        if (! exists $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM_non_multi'} ) {
            $data_ref->{$imp_gene}->{$expr_type2}->{'RPKM_non_multi'} = 0.00;
        }
    }

    # Get concise strings for individual RPKM data, both simple and non-multi.:
    my $L3_rpkm    = $data_ref->{$imp_gene}->{'L3'}->{'RPKM'};
    my $L4_rpkm    = $data_ref->{$imp_gene}->{'L4'}->{'RPKM'};
    my $nhr67_rpkm = $data_ref->{$imp_gene}->{'nhr'}->{'RPKM'};
    my $N2_rpkm    = $data_ref->{$imp_gene}->{'N2'}->{'RPKM'};

    my $L3_rpkm_nm    = $data_ref->{$imp_gene}->{'L3'}->{'RPKM_non_multi'};
    my $L4_rpkm_nm    = $data_ref->{$imp_gene}->{'L4'}->{'RPKM_non_multi'};     
    my $nhr67_rpkm_nm = $data_ref->{$imp_gene}->{'nhr'}->{'RPKM_non_multi'};
    my $N2_rpkm_nm    = $data_ref->{$imp_gene}->{'N2'}->{'RPKM_non_multi'};

    # Get a concise string variable for maximum RPKM of whatever values were considered important.
    # The default is to get the maximum of wild-type L3 and L4 (because this is the default "important").
    # But making this respond to the actual choice of important values means that the table always ends up 
    #     sorting by the maximum value of what was being asked about.
    # Also get a maximum for non-multi. RPKM.
    # Conceivably, this could come from a different source than the ordinary maximum, but in practice that is unlikely.

    my @important_rpkms = ();
    foreach my $imp_dat (sort @important_data) {
        push @important_rpkms, $data_ref->{$imp_gene}->{$imp_dat}->{'RPKM'};
    }
    my @important_rpkm_nms = ();
    foreach my $imp_dat (sort @important_data) {
        push @important_rpkm_nms, $data_ref->{$imp_gene}->{$imp_dat}->{'RPKM_non_multi'};
    }

    my $max_imp_rpkm = max( @important_rpkms );
    $data_ref->{$imp_gene}->{'max_imp_RPKM'} = $max_imp_rpkm;
    my $max_imp_rpkm_nm = max( @important_rpkm_nms );
    $data_ref->{$imp_gene}->{'max_imp_RPKM_nm'} = $max_imp_rpkm_nm;

    # However, we always also want L3/L4 max.,  for ratios with N2, nhr-67(-), etc.  So:
    my $max_wtLC_rpkm = max ( $data_ref->{$imp_gene}->{'L3'}->{'RPKM'}, $data_ref->{$imp_gene}->{'L4'}->{'RPKM'} );
    $data_ref->{$imp_gene}->{'max_wtLC_RPKM'} = $max_wtLC_rpkm;
    my $max_wtLC_rpkm_nm = max( $data_ref->{$imp_gene}->{'L3'}->{'RPKM_non_multi'}, $data_ref->{$imp_gene}->{'L4'}->{'RPKM_non_multi'} );
    $data_ref->{$imp_gene}->{'max_wtLC_RPKM_nm'} = $max_wtLC_rpkm_nm;

    # Get log[10] versions of these worked up and stored as text-lineseg.
    # Can use safely_divide_RPKMs, since getting the log of X is equivalent to getting the log of X/1!
    #     In other words, the job I might want to do with safely_take_log is a subset of what safely_divide_RPKMs does.
    #     So, just invoke that with a simple hack...

    my $L3_log10          = safely_divide_RPKMs($L3_rpkm,          1, 0.01, 'log', 10);
    my $L4_log10          = safely_divide_RPKMs($L4_rpkm,          1, 0.01, 'log', 10);
    my $nhr_log10         = safely_divide_RPKMs($nhr67_rpkm,       1, 0.01, 'log', 10);
    my $N2_log10          = safely_divide_RPKMs($N2_rpkm,          1, 0.01, 'log', 10);

    my $L3_nm_log10       = safely_divide_RPKMs($L3_rpkm_nm,       1, 0.01, 'log', 10);
    my $L4_nm_log10       = safely_divide_RPKMs($L4_rpkm_nm,       1, 0.01, 'log', 10);
    my $nhr_nm_log10      = safely_divide_RPKMs($nhr67_rpkm_nm,    1, 0.01, 'log', 10);
    my $N2_nm_log10       = safely_divide_RPKMs($N2_rpkm_nm,       1, 0.01, 'log', 10);

    my $max_imp_log10     = safely_divide_RPKMs($max_imp_rpkm,     1, 0.01, 'log', 10);
    my $max_imp_nm_log10  = safely_divide_RPKMs($max_imp_rpkm_nm,  1, 0.01, 'log', 10);

    my $max_wtLC_log10    = safely_divide_RPKMs($max_wtLC_rpkm,    1, 0.01, 'log', 10);
    my $max_wtLC_nm_log10 = safely_divide_RPKMs($max_wtLC_rpkm_nm, 1, 0.01, 'log', 10);

    # Get the [max. imp.]/[N2] ratio data, in human-readable plain and log[10] form:
    my $max_imp_N2_ratio           = safely_divide_RPKMs($max_imp_rpkm, $N2_rpkm, 0.01);
    my $max_imp_N2_ratio_log10     = safely_divide_RPKMs($max_imp_rpkm, $N2_rpkm, 0.01, 'log', 10);
    
    # Get the [max. imp. (non-multiread)]/[N2] ratio data, in human-readable plain and log[10] form:
    my $max_imp_nm_N2_ratio           = safely_divide_RPKMs($max_imp_rpkm_nm, $N2_rpkm, 0.01);
    my $max_imp_nm_N2_ratio_log10     = safely_divide_RPKMs($max_imp_rpkm_nm, $N2_rpkm, 0.01, 'log', 10);

    # Get the [max. wt LC]/[N2] ratio data, in human-readable plain and log[10] form:
    my $max_wtLC_N2_ratio          = safely_divide_RPKMs($max_wtLC_rpkm, $N2_rpkm, 0.01);
    my $max_wtLC_N2_ratio_log10    = safely_divide_RPKMs($max_wtLC_rpkm, $N2_rpkm, 0.01, 'log', 10);

    # Get the [max wt LC (non-multiread)]/[N2] ratio data, in human-readable plain and log[10] form:
    my $max_wtLC_nm_N2_ratio       = safely_divide_RPKMs($max_wtLC_rpkm_nm, $N2_rpkm, 0.01); 
    my $max_wtLC_nm_N2_ratio_log10 = safely_divide_RPKMs($max_wtLC_rpkm_nm, $N2_rpkm, 0.01, 'log', 10); 

    # Get the [wt L4]/[wt L3] ratio data, with simple and non-multi. RPKMs:
    my $L4_L3_ratio          = safely_divide_RPKMs($L4_rpkm,    $L3_rpkm,    0.01);
    my $L4_L3_ratio_log10    = safely_divide_RPKMs($L4_rpkm,    $L3_rpkm,    0.01, 'log', 10);
    my $L4_L3_nm_ratio       = safely_divide_RPKMs($L4_rpkm_nm, $L3_rpkm_nm, 0.01);
    my $L4_L3_nm_ratio_log10 = safely_divide_RPKMs($L4_rpkm_nm, $L3_rpkm_nm, 0.01, 'log', 10);

    # Get the [nhr-67(-) L4]/[wt L3] ratio data, with simple and non-multi. RPKMs:
    my $nhr67_L3_ratio          = safely_divide_RPKMs($nhr67_rpkm,    $L3_rpkm,    0.01);
    my $nhr67_L3_ratio_log10    = safely_divide_RPKMs($nhr67_rpkm,    $L3_rpkm,    0.01, 'log', 10);
    my $nhr67_L3_nm_ratio       = safely_divide_RPKMs($nhr67_rpkm_nm, $L3_rpkm_nm, 0.01);
    my $nhr67_L3_nm_ratio_log10 = safely_divide_RPKMs($nhr67_rpkm_nm, $L3_rpkm_nm, 0.01, 'log', 10);

    # Get the [nhr-67(-) L4]/[wt L4] ratio data, with simple and non-multi. RPKMs:
    my $nhr67_L4_ratio          = safely_divide_RPKMs($nhr67_rpkm,    $L4_rpkm,    0.01);
    my $nhr67_L4_ratio_log10    = safely_divide_RPKMs($nhr67_rpkm,    $L4_rpkm,    0.01, 'log', 10);
    my $nhr67_L4_nm_ratio       = safely_divide_RPKMs($nhr67_rpkm_nm, $L4_rpkm_nm, 0.01);
    my $nhr67_L4_nm_ratio_log10 = safely_divide_RPKMs($nhr67_rpkm_nm, $L4_rpkm_nm, 0.01, 'log', 10);

    # Store these data individually, so that I can later use them to sort the output:
    $data_ref->{'final'}->{$imp_gene}->{'max_imp_vs_N2_ratio'}       = $max_imp_N2_ratio;
    $data_ref->{'final'}->{$imp_gene}->{'max_imp_vs_N2_ratio_log10'} = $max_imp_N2_ratio_log10;

    # After all scalars computed, then round them off (by default, but optionally not).
    # Note that this rounding does *not* touch values stored in $data_ref, 
    #     such as: $data_ref->{$imp_gene}->{'max_wtLC_RPKM_nm'} = $max_wtLC_rpkm_nm;
    # But it will affect the following scalars, which will then be marshalled into arrays for printing.

    if (! $no_rounding) {
        # Basic data:
        # No need to round $L3_rpkm, $L4_rpkm, $nhr67_rpkm, or $N2_rpkm -- they're prerounded to 0.01!

        $L3_rpkm_nm      = sprintf "%.2f", $L3_rpkm_nm;
        $L4_rpkm_nm      = sprintf "%.2f", $L4_rpkm_nm;
        $nhr67_rpkm_nm   = sprintf "%.2f", $nhr67_rpkm_nm;
        $N2_rpkm_nm      = sprintf "%.2f", $N2_rpkm_nm;

        # Log(10) of basic data:
        $L3_log10        = sprintf "%.3f", $L3_log10;
        $L4_log10        = sprintf "%.3f", $L4_log10;
        $nhr_log10       = sprintf "%.3f", $nhr_log10;
        $N2_log10        = sprintf "%.3f", $N2_log10;

        $L3_nm_log10     = sprintf "%.3f", $L3_nm_log10;
        $L4_nm_log10     = sprintf "%.3f", $L4_nm_log10;
        $nhr_nm_log10    = sprintf "%.3f", $nhr_nm_log10;
        $N2_nm_log10     = sprintf "%.3f", $N2_nm_log10;

        # Maxima:
        $max_imp_rpkm     = sprintf "%.2f", $max_imp_rpkm;   
        $max_imp_rpkm_nm  = sprintf "%.2f", $max_imp_rpkm_nm;
        $max_wtLC_rpkm    = sprintf "%.2f", $max_wtLC_rpkm;
        $max_wtLC_rpkm_nm = sprintf "%.2f", $max_wtLC_rpkm_nm;

        # Log(10) of maxima:
        $max_imp_log10     = sprintf "%.3f", $max_imp_log10;
        $max_imp_nm_log10  = sprintf "%.3f", $max_imp_nm_log10;
        $max_wtLC_log10    = sprintf "%.3f", $max_wtLC_log10;
        $max_wtLC_nm_log10 = sprintf "%.3f", $max_wtLC_nm_log10;

        # [Maximum 'important']/[N2]:
        $max_imp_N2_ratio           = sprintf "%.2f", $max_imp_N2_ratio;
        $max_imp_N2_ratio_log10     = sprintf "%.3f", $max_imp_N2_ratio_log10;
        $max_imp_nm_N2_ratio        = sprintf "%.2f", $max_imp_nm_N2_ratio;
        $max_imp_nm_N2_ratio_log10  = sprintf "%.3f", $max_imp_nm_N2_ratio_log10;

        # [Max. LC]/[N2]:
        $max_wtLC_N2_ratio          = sprintf "%.2f", $max_wtLC_N2_ratio;
        $max_wtLC_N2_ratio_log10    = sprintf "%.3f", $max_wtLC_N2_ratio_log10;
        $max_wtLC_nm_N2_ratio       = sprintf "%.2f", $max_wtLC_nm_N2_ratio;
        $max_wtLC_nm_N2_ratio_log10 = sprintf "%.3f", $max_wtLC_nm_N2_ratio_log10;

        # Other ratios:
        $L4_L3_ratio             = sprintf "%.2f", $L4_L3_ratio;
        $L4_L3_ratio_log10       = sprintf "%.3f", $L4_L3_ratio_log10;
        $L4_L3_nm_ratio          = sprintf "%.2f", $L4_L3_nm_ratio;
        $L4_L3_nm_ratio_log10    = sprintf "%.3f", $L4_L3_nm_ratio_log10;

        $nhr67_L3_ratio          = sprintf "%.2f", $nhr67_L3_ratio;
        $nhr67_L3_ratio_log10    = sprintf "%.3f", $nhr67_L3_ratio_log10;
        $nhr67_L3_nm_ratio       = sprintf "%.2f", $nhr67_L3_nm_ratio;
        $nhr67_L3_nm_ratio_log10 = sprintf "%.3f", $nhr67_L3_nm_ratio_log10;

        $nhr67_L4_ratio          = sprintf "%.2f", $nhr67_L4_ratio;
        $nhr67_L4_ratio_log10    = sprintf "%.3f", $nhr67_L4_ratio_log10;
        $nhr67_L4_nm_ratio       = sprintf "%.2f", $nhr67_L4_nm_ratio;
        $nhr67_L4_nm_ratio_log10 = sprintf "%.3f", $nhr67_L4_nm_ratio_log10;
    }

    # Archive *all* column header texts, and *all* data wrangled above, as arrayrefs.
    # This makes both output and rearranging the columns simple.

    my @header_columns = qw( Gene

                             L3            L4             nhr-67            N2            max_wtLC            max_wtLC/N2
                             L3_nm         L4_nm          nhr-67_nm         N2_nm         max_wtLC_nm         max_wtLC_nm/N2_nm
                             log10(L3)     log10(L4)      log10(nhr-67)     log10(N2)     log10(max_wtLC)     log10(max_wtLC/N2)
                             log10(L3_nm)  log10(L4_nm)   log10(nhr-67_nm)  log10(N2_nm)  log10(max_wtLC_nm)  log10(max_wtLC_nm/N2)

                             L4/L3         log10(L4/L3)   L4_nm/L3_nm       log10(L4_nm/L3_nm)
                             nhr/L3        log10(nhr/L3)  nhr_nm/L3_nm      log10(nhr_nm/L3_nm)
                             nhr/L4        log10(nhr/L4)  nhr_nm/L4_nm      log10(nhr_nm/L4_nm)
    );

    my @auxillary_headers = qw( maxImp     maxImp_nm         log10(maxImp)  log10(maxImp_nm)
                                maxImp/N2  log10(maxImp/N2)  maxImp_nm/N2   log10(maxImp_nm/N2)  );
    if ($auxillary) { 
        push @header_columns, @auxillary_headers;
    }

    my $header_columns_ref = \@header_columns;
    $data_ref->{'header_columns'} = $header_columns_ref;

    my @data_to_print = ( $L3_rpkm,     $L4_rpkm,     $nhr67_rpkm,    $N2_rpkm,     $max_wtLC_rpkm,     $max_wtLC_N2_ratio,
                          $L3_rpkm_nm,  $L4_rpkm_nm,  $nhr67_rpkm_nm, $N2_rpkm_nm,  $max_wtLC_rpkm_nm,  $max_wtLC_nm_N2_ratio, 
                          $L3_log10,    $L4_log10,    $nhr_log10,     $N2_log10,    $max_wtLC_log10,    $max_wtLC_N2_ratio_log10, 
                          $L3_nm_log10, $L4_nm_log10, $nhr_nm_log10,  $N2_nm_log10, $max_wtLC_nm_log10, $max_wtLC_nm_N2_ratio_log10, 

                          $L4_L3_ratio,       $L4_L3_ratio_log10,    $L4_L3_nm_ratio,    $L4_L3_nm_ratio_log10,
                          $nhr67_L3_ratio,    $nhr67_L3_ratio_log10, $nhr67_L3_nm_ratio, $nhr67_L3_nm_ratio_log10,
                          $nhr67_L4_ratio,    $nhr67_L4_ratio_log10, $nhr67_L4_nm_ratio, $nhr67_L4_nm_ratio_log10,      );

    my @auxillary_data = ( $max_imp_rpkm,     $max_imp_rpkm_nm,        $max_imp_log10,       $max_imp_nm_log10, 
                           $max_imp_N2_ratio, $max_imp_N2_ratio_log10, $max_imp_nm_N2_ratio, $max_imp_nm_N2_ratio_log10 );

    if ($auxillary) { 
        push @data_to_print, @auxillary_data;
    }
    my $data_to_print_ref = \@data_to_print;
    $data_ref->{'final'}->{$imp_gene}->{'data_to_print'} = $data_to_print_ref;
}

my @final_genes = sort { $data_ref->{'final'}->{$b}->{'max_imp_vs_N2_ratio'} 
                         <=> $data_ref->{'final'}->{$a}->{'max_imp_vs_N2_ratio'} } 
                  sort { $data_ref->{'final'}->{$b}->{'max_imp_vs_N2_ratio_log10'}
                         <=> $data_ref->{'final'}->{$a}->{'max_imp_vs_N2_ratio_log10'} }
                  sort { $data_ref->{$b}->{'max_imp_RPKM'} 
                         <=> $data_ref->{$a}->{'max_imp_RPKM'} } 
                  keys %{ $data_ref->{'final'} };


if ($header_text) {
    my $header = join "\t", @{ $data_ref->{'header_columns'} };
    print "$header\n";
}

foreach my $final_gene (@final_genes) { 
    print $final_gene,
          "\t",
          (join "\t", @{ $data_ref->{'final'}->{$final_gene}->{'data_to_print'} } ),
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

