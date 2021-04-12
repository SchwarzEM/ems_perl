#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $infile        = q{};
my $sign          = q{};
my $log2fc_thresh = q{};
my $fdr_thresh    = q{};

if (! $ARGV[3]) {
    die "Format: filter_immunoregs_12apr2021.pl [infile] [sign (up or down)] [log2FC_threshold] [FDR threshold] > [output]\n"; 
}

$infile        = $ARGV[0] if $ARGV[0] ;
$sign          = $ARGV[1] if $ARGV[1] ;
$log2fc_thresh = $ARGV[2] if $ARGV[2] ;
$fdr_thresh    = $ARGV[3] if $ARGV[3] ;

if ( ( $sign ne 'up' ) and ( $sign ne 'down' ) ) {
    die "Sign must be either 'up' or 'down', but is: $sign\n";
}

if (! looks_like_number($log2fc_thresh) ) {
    die "log2FC_threshold is not a number: $log2fc_thresh\n";
}

if ( (! looks_like_number($fdr_thresh) ) or ( $fdr_thresh < 0 ) or ( $fdr_thresh > 1 ) ) {
    die "FDR_threshold needs to be a number between 0 and 1, but is: $fdr_thresh\n";
}

if ( ( $log2fc_thresh < 0 ) and ( $sign eq 'up' ) ) {
    warn "Note that log2FC_threshold is negative ($log2fc_thresh) but sign is $sign\n";
}

if ( ( $log2fc_thresh > 0 ) and ( $sign eq 'down' ) ) {
    warn "Note that log2FC_threshold is positive ($log2fc_thresh) but sign is $sign\n";
}

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A Gene \t/xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A (?: [^\t]* \t){45} ([^\t]*) \t ([^\t]*) \t ([^\t]*) \t ([^\t]*) \t .* \z/xms ) {
        my $int_logFC     = $1;
        my $int_FDR       = $2;

        my $non_int_logFC = $3;
        my $non_int_FDR   = $4;

        if (     ( $int_logFC =~ /\S/ ) 
             and ( $int_FDR =~ /\S/ ) 
             and ( $int_FDR <= $fdr_thresh ) 
             and (    ( ( $int_logFC >= $log2fc_thresh ) and ($sign eq 'up' )   )
                   or ( ( $int_logFC <= $log2fc_thresh ) and ($sign eq 'down' ) )
                 )
            ) {
            print "$input\n";
        }
        elsif (     ( $non_int_logFC =~ /\S/ )
                and ( $non_int_FDR =~ /\S/ ) 
                and ( $non_int_FDR <= $fdr_thresh ) 
                and (    ( ( $non_int_logFC >= $log2fc_thresh ) and ($sign eq 'up' )   )
                      or ( ( $non_int_logFC <= $log2fc_thresh ) and ($sign eq 'down' ) )
                    )
              ) {
            print "$input\n";
        }
    }
    else {
        die "Cannot parse: $input\n";
    }
}

close $INFILE;

