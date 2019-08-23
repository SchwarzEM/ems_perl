#!/usr/bin/env perl

# coverage_fastq_residues.pl -- Erich Schwarz <ems394@cornell.edu>, 8/23/2019. 
# Purpose: get statistics from a FASTQ files or stream (max., min., mean, median, N50 etc., optional coverages).

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);

my $infile    = q{};
my $nt_count  = 0;
my @sizes     = ();
my @n_ranks   = (10, 20, 30, 40, 50, 60, 70, 80, 90);
my $genome    = q{};
my @coverages = ();
my $i         = 0;

my $help;

GetOptions ( 'infile=s'       => \$infile,
             'genome=s'       => \$genome,
             'coverages=s{,}' => \@coverages,
             'help'           => \$help,   );

if (    $help 
     or (! $infile ) 
     or ( $genome     and (! @coverages) )
     or ( (! $genome) and @coverages     )
   ) { 
    die "Format: count_simple_fastq_residues.pl\n",
        "    --infile|-i     <input stream/file>\n",
        "    --genome|-g     <genome size for coverage; requires 1+ coverages; values of xG, xM, or xK parse as Gb, Mb, or kb>\n",
        "    --coverages|-c  <1+ coverages to compute; requires a genome size>\n",
        "    --help|-h       [print this message]\n",
        ;
}

foreach my $coverage (@coverages) {
    if (! looks_like_number($coverage) ) {
        die "Non-numeric coverage value: $coverage\n";
    }
    if ( $coverage <= 0 ) {
        die "Non-positive coverage value: $coverage\n";
    }
}   

if ( $genome and ( $genome =~ /\A (\S+) (G|g|M|m|K|k) \z/xms ) ) {
    my $value  = $1;
    my $suffix = $2;

    if (! looks_like_number($value) ) {
        die "Cannot parse numerical value in genome size: $genome\n";
    }

    if ( $suffix =~ /\A G|g \z/xms ) {
        $genome = $value * 1_000_000_000;
    }
    elsif ( $suffix =~ /\A M|m \z/xms ) {
       	$genome = $value * 1_000_000;
    }
    elsif ( $suffix =~ /\A K|k \z/xms ) {
        $genome = $value * 1_000;
    }
    else {
        die "Cannot parse suffix in genome size: $genome\n";
    }

    if (! looks_like_number($genome) ) {
       	die "Even after parsing, non-numerical genome size: $genome\n";
    }
}

if ( $genome and (! looks_like_number($genome) ) ) {
    die "Non-numerical genome size: $genome\n";
}

if ( $genome and ( $genome <= 0 ) ) {
    die "Non-positive genome size: $genome\n";
}

my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}
while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    $i++;
    # Use !(modulus) function to get sequence from every 2cd line of the 4-line records.
    if ( ( ( $i % 4 ) == 2 ) and ( $input =~ / \A \s* [ACGTNacgt]+ \s* \z /xms ) ) { 
        $input =~ s/\s//g;
        $nt_count += length($input);
        push @sizes, $nt_count;
        $nt_count = 0;
    }
    if ( ( ( $i % 4 ) == 2 ) and ( $input !~ / \A \s* [ACGTNacgt]+ \s* \z /xms ) ) {
        die "Unexpected sequence: $input\n";
    }
}
close $INPUT_FILE;

# Sort in descending numerical order, so that homebrewed N_value (N50, etc.) subroutine can work.
@sizes = sort { $b <=> $a } @sizes;

# Sort in descending numerical order for optional coverage analysis.
@coverages = sort { $b <=> $a } @coverages ;

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@sizes);

my $sum     = $stat->sum();    # Total nt of sequence.
my $rawsum = $sum;
$sum = commify($sum);

my $reads = $stat->count();  # Total number of reads.
$reads = commify($reads);

my $mean    = $stat->mean();
$mean       = sprintf("%.2f", $mean);

my $std_dev = $stat->standard_deviation();
$std_dev    = sprintf("%.2f", $std_dev);

my $min     = $stat->min();
$min        = commify($min);

my $max     = $stat->max();
$max        = commify($max);

my $median = $stat->median();
$median    = sprintf("%.2f", $median);

print "Reads:   $reads\n\n";
print "Total:   $sum nt\n\n";
print "Mean:    $mean nt\n";
print "Median:  $median nt\n";
print "Max:     $max nt\n";
print "Min:     $min nt\n\n";
foreach my $n_rank (@n_ranks) {
    my $n_value = &get_n_value(\@sizes,$rawsum,$n_rank);

    # Get a reasonable decimal fraction (two digits).
    $n_value = sprintf("%.2f", $n_value );

    # Round to nearest digit (in Perl, this means adding 0.5, then 'int'.
    $n_value = $n_value + 0.5;
    $n_value = int($n_value);

    # Make human-readable, and print.
    $n_value = commify($n_value);
    print "N", "$n_rank:     $n_value nt\n";
}

if ( $genome and @coverages ) {
    my $gen_read = commify($genome);
    my $cov_header = "\nGenome size: $gen_read nt\n\n";
    foreach my $coverage (@coverages) {
        my $cov_value = &get_cov_value(\@sizes,$genome,$coverage);

        # Get a reasonable decimal fraction (two digits).
        $cov_value = sprintf("%.2f", $cov_value );

        # Round to nearest digit (in Perl, this means adding 0.5, then 'int'.
        $cov_value = $cov_value + 0.5;
        $cov_value = int($cov_value);

        # Make human-readable.
        $cov_value = commify($cov_value);

        # Have one space before all coverages.
        print $cov_header if $cov_header;
        $cov_header = q{};

        # Print human-readable coverage.
        print "$coverage", "x coverage: $cov_value nt\n";
    }
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

# Assumes array sorted in ascending values, with precomputed total value:
sub get_n_value { 
    my ($_sizes_ref, $_sum, $_n_rank) = @_;
    my $_n_sum = ($_sum * ( $_n_rank / 100));
    my $_tally = 0;
    my @_sizes = @{ $_sizes_ref };
    foreach my $_value (@_sizes) { 
        $_tally += $_value;
        if ($_tally > $_n_sum) { 
            return $_value;
        }
    }
    return;
}

sub get_cov_value {
    my ($_sizes_ref, $_genome, $_coverage) = @_;
    my $_cov_size = ($_genome * $_coverage);
    my $_tally = 0;
    my @_sizes = @{ $_sizes_ref };
    foreach my $_value (@_sizes) {
        $_tally += $_value;
        if ($_tally > $_cov_size) {
            return $_value;
        }
    }
    return;
}
