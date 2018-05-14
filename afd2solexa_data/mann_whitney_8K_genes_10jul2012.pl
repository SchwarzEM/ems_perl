#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Statistics::Test::WilcoxonRankSum;

my $infile           = q{};
my $pattern          = q{};
my @marked_genes     = ();
my @non_marked_genes = ();
my $help;

GetOptions ( 'infile=s'  => \$infile,
             'pattern=s' => \$pattern,
             'help'      => \$help, ); 

if ($help or (! $infile) ) {
    print "Format: mann_whitney_8K_genes_10jul2012.pl\n",
          "            --infile|-i    [input text file or '-']\n",
          "            --pattern|-p   [a pattern which will be used to select marked genes in the total set]\n",
          "            --help|-h      [print this message]\n",
          ;
    exit;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}

# Keep quotable version of pattern:
my $text_pattern = $pattern;

# Do special Perl quoting to make the pattern work better in a regex:
$pattern = qr/$pattern/;

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S+ \t (\S+) \t (\S+) \t (\S*) /xms ) { 
        my $rpkm_l3    = $1;
        my $rpkm_nonl3 = $2;
        my $marker     = $3;
        if ( $rpkm_l3 == 0 ) { 
            $rpkm_l3 = 0.01;
        }
        if ( $rpkm_nonl3 == 0 ) {
            $rpkm_nonl3 = 0.01;
        }
        my $ratio_nonl3_l3 = ( $rpkm_nonl3 / $rpkm_l3 );
        # Omit '/xms' to make pattern matches more reliable.
        if ( $marker =~ /$pattern/ ) { 
            push @marked_genes, $ratio_nonl3_l3;
        }
        else { 
            push @non_marked_genes, $ratio_nonl3_l3;
        }
    }
    else { 
       die "Can't parse input: $input\n";
    }
}
close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();

$wilcox_test->load_data(\@marked_genes, \@non_marked_genes);

print "Dataset 1: genes marked with pattern \"$text_pattern\"\n";
print "Dataset 2: other genes\n";

my $prob = $wilcox_test->probability();
my $pf    = sprintf '%f', $prob;
print "$pf\n";

my $pstatus = $wilcox_test->probability_status();
print "$pstatus\n";

my $psum = $wilcox_test->summary();
print "$psum\n";

