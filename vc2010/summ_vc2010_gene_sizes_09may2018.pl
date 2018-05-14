#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Statistics::Descriptive;

my @all_sizes_array   = ();
my @prot_sizes_array  = ();
my @ncrna_sizes_array = ();

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t gene \t (\d+) \t (\d+) \t /xms ) { 
        my $start_nt = $1;
        my $end_nt   = $2;
        if ( $end_nt < $start_nt ) {
            die "Cannot parse gene coordinates: $input\n";
        }

        my $size = ( $end_nt - $start_nt + 1);

        push @all_sizes_array, $size;
        if ( $input =~ / protein_coding /xms ) {
            push @prot_sizes_array, $size;
        }
        else {
            push @ncrna_sizes_array, $size;
        }
    }
    else { 
        die "Cannot parse: $input\n"
    }
}

my $all_sizes = Statistics::Descriptive::Full->new();
$all_sizes->add_data(@all_sizes_array);

my $all_sizes_count  = $all_sizes->count();  # Total number of genes
$all_sizes_count     = commify($all_sizes_count);

my $all_sizes_mean   = $all_sizes->mean();
$all_sizes_mean      = sprintf("%.1f", $all_sizes_mean);
$all_sizes_mean      = commify($all_sizes_mean);

my $all_sizes_median = $all_sizes->median();
$all_sizes_median    = sprintf("%.1f", $all_sizes_median);
$all_sizes_median    = commify($all_sizes_median);

my $all_sizes_max    = $all_sizes->max();
$all_sizes_max       = commify($all_sizes_max);

my $all_sizes_min    = $all_sizes->min();
$all_sizes_min       = commify($all_sizes_min);

my $prot_sizes = Statistics::Descriptive::Full->new();
$prot_sizes->add_data(@prot_sizes_array);

my $prot_sizes_count  = $prot_sizes->count();  # Total number of genes
$prot_sizes_count     = commify($prot_sizes_count);
        
my $prot_sizes_mean   = $prot_sizes->mean();
$prot_sizes_mean      = sprintf("%.1f", $prot_sizes_mean);
$prot_sizes_mean      = commify($prot_sizes_mean);
        
my $prot_sizes_median = $prot_sizes->median();
$prot_sizes_median    = sprintf("%.1f", $prot_sizes_median);
$prot_sizes_median    = commify($prot_sizes_median);

my $prot_sizes_max    = $prot_sizes->max();
$prot_sizes_max       = commify($prot_sizes_max);

my $prot_sizes_min    = $prot_sizes->min();
$prot_sizes_min       = commify($prot_sizes_min);

my $ncrna_sizes = Statistics::Descriptive::Full->new();
$ncrna_sizes->add_data(@ncrna_sizes_array);

my $ncrna_sizes_count  = $ncrna_sizes->count();  # Total number of genes
$ncrna_sizes_count     = commify($ncrna_sizes_count);
    
my $ncrna_sizes_mean   = $ncrna_sizes->mean();
$ncrna_sizes_mean      = sprintf("%.1f", $ncrna_sizes_mean);
$ncrna_sizes_mean      = commify($ncrna_sizes_mean);
    
my $ncrna_sizes_median = $ncrna_sizes->median();
$ncrna_sizes_median    = sprintf("%.1f", $ncrna_sizes_median);
$ncrna_sizes_median    = commify($ncrna_sizes_median);

my $ncrna_sizes_max    = $ncrna_sizes->max();
$ncrna_sizes_max       = commify($ncrna_sizes_max);

my $ncrna_sizes_min    = $ncrna_sizes->min();
$ncrna_sizes_min       = commify($ncrna_sizes_min);

print "All genes --\n\n";
print "Count:            $all_sizes_count\n";
print "Mean size, nt:    $all_sizes_mean\n";
print "Median size, nt:  $all_sizes_median\n";
print "Max size, nt:     $all_sizes_max\n";
print "Min size, nt:     $all_sizes_min\n";
print "\n";
print "Protein-coding genes --\n\n";
print "Count:            $prot_sizes_count\n";   
print "Mean size, nt:    $prot_sizes_mean\n";      
print "Median size, nt:  $prot_sizes_median\n";
print "Max size, nt:     $prot_sizes_max\n";       
print "Min size, nt:     $prot_sizes_min\n";       
print "\n";
print "Non-protein-coding genes --\n\n";
print "Count:            $ncrna_sizes_count\n";
print "Mean size, nt:    $ncrna_sizes_mean\n";
print "Median size, nt:  $ncrna_sizes_median\n";
print "Max size, nt:     $ncrna_sizes_max\n";
print "Min size, nt:     $ncrna_sizes_min\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


