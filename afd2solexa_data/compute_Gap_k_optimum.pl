#!/usr/bin/env perl

# compute_Gap_k_optimum.pl -- Erich Schwarz <ems394@cornell.edu>, 3/2/2013.
# Purpose: given a 'vec' or 'sim'-compatible data matrix, compute its Gap optimum.

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Statistics::Gap;

my $output_prefix = q{};
my $type          = 'vec';
my $input_matrix  = q{};
my $clust_meth    = 'rbr';
my $criterion     = 'i2';
my $upper_bound   = 0;
my $number_rep    = 0;
my $rep_opt       = 'rep';
my $perc_conf     = 90;
my $precision     = 4;
my $seed          = q{};
my $help;

my %ok_clust_meths = ( 'rb' => 1, 
                       'rbr' => 1, 
                       'direct' => 1,
                       'agglo' => 1,
                       'bagglo' => 1, );

my %ok_criteria = ( 'i1' => 1, 
                    'i2' => 1,
                    'e1' => 1,
                    'h1' => 1,
                    'h2' => 1, );


GetOptions ( 'output_prefix=s' => \$output_prefix,
             'type=s'          => \$type,
             'input_matrix=s'  => \$input_matrix,
             'meth_clust=s'    => \$clust_meth,
             'criterion=s'     => \$criterion,
             'upper_bound=i'   => \$upper_bound,
             'number_rep=i'    => \$number_rep,
             'rep_opt=s'       => \$rep_opt,
             'perc_conf=f'     => \$perc_conf,
             'val_precision=i' => \$precision,
             'seed=i'          => \$seed,
             'help'            => \$help, );

if ($help 
       or ( ( $type ne 'vec' ) and ( $type ne 'sim' ) )  
       or (! -r $input_matrix ) 
       or (! $ok_clust_meths{$clust_meth} ) or ( ( $clust_meth eq 'bagglo' ) and ( $type ne 'vec' ) ) 
       or (! $ok_criteria{$criterion} )
       or ( ( $upper_bound != int $upper_bound ) or ( $upper_bound < 1 ) ) 
       or ( ( $number_rep  != int $number_rep  ) or ( $number_rep  < 1 ) ) 
       or ( ( $rep_opt ne 'rep' ) and ( $rep_opt ne 'ref' ) ) 
       or ( ( $perc_conf <= 0 ) or ( $perc_conf > 100 ) ) 
       or ( ( $precision != int $precision ) or ( $precision < 1 ) ) 
    ) { 
    die "Format: compute_Gap_k_optimum.pl\n",
        "        --output_prefix|-o   [prefix for naming intermediate files and .dat files]\n",
        "        --type|-t            [type of space in which clustering will be performed: 'vec' (vector) or 'sim'; default 'vec']\n",
        "        --input_matrix|-i    [path to input matrix]\n",
        "                                 [N.B.: for more details on following options, see 'perldoc Statistics::Gap' and refs. therein]\n",
        "        --meth_clust|-m      [clustering method: 'rb', 'rbr', 'direct', 'agglo', 'bagglo' (if 'vec'); default is 'rbr']\n",
        "        --criterion|-c       [criterion function for finding clustering solutions: 'i1', 'i2', 'e1', 'h1', 'h2'; default is 'i2']\n",
        "        --upper_bound|-u     [user-specified upper bound for allowable number of clusters]\n",
        "        --number_rep|-n      [number of replicates/references to be generated]\n",
        "        --rep_opt|-r         [generate R replicates from a reference ('rep') or generate B references ('ref'); default 'rep']\n",
        "        --perc_conf|-p       [percentage confidence to be reported in the log file, for parametric bootstrapping: default is 90]\n",
        "        --val_precision|-v   [precision to be used while generating the reference distribution: default is 4]\n",
        "        --seed|-s            [seed to be used with the random number generator (optional argument)]\n",
        ;
}

my $basename = basename($input_matrix);

my $predictedk = &gap($output_prefix,
           $type,
           $input_matrix,
           $clust_meth,
           $criterion,
           $upper_bound,
           $number_rep,
           $rep_opt,
           $perc_conf,
           $precision,
           $seed,
); 

print "\n",
      "For data in $input_matrix, with:",
      "    clustering type    $type\n",
      "    clustering method  $clust_meth\n",
      "    criterion function $criterion\n",
      "    upper bound for k  $upper_bound\n",
      "    replicates         $number_rep\n",
      "    rep. option        $rep_opt\n",
      "    percent. conf.     $perc_conf\n",
      "    and precision      $precision\n",
      "\n",
      "    predicted optimal k value is: $predictedk\n",
      "\n",
      ;

