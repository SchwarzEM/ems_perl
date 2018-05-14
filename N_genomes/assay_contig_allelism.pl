#!/usr/bin/env perl

# assay_contig_allelism.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/3/2009.
# Purpose: list the effective nt. count, with mode / median / mean per-nt-copy no., for Shiaw-Pyng Yang's contigs.

# N.B.: an earlier version that counted all values and did mean, median, etc. with
#    Statistics::Descriptive totally crashed system memory on genome-size files.

use strict;
use warnings;
use Getopt::Long;

my $tabs;
my $fixed;
my @output_vals = ();

GetOptions ( "tabs"  => \$tabs,
             "fixed" => \$fixed, ); 

if ((! $tabs ) and (! $fixed)) { 
    die "Format: assay_contig_allelism.pl --tabs|--fixed (or -t|-f) [files/STDIN]\n";
}
if ($tabs and $fixed) {
    die "Either tab-delimited or fixed-width output, but not both!\n";
}

my $contig = q{};
my $contig_info_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\w+) \. /xms ) { 
        $contig = $1;
    } 
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse header line: $input\n";
    }
    elsif ( $input =~ /\A \s* \d+ (?: \s+ \d+)* \s* \z /xms ) { 
        $input =~ s/\A\s+//;
        $input =~ s/\s+\z//;
        my @numbers = split /\s+/, $input;
        foreach my $i (@numbers) { 
            $contig_info_ref->{$contig}->{'total'} += 1;
            $contig_info_ref->{$contig}->{'values'}->{$i} += 1;
        }
    }
    elsif ( $input =~ / \S /xms ) { 
        die "Can't parse information line: $input\n";
    }
}

@output_vals = qw( Contig  Nt(*)   Mode 
                   Mean    Frac_1  Frac_2 );

print_output_line(\@output_vals, $tabs, $fixed);

foreach my $id (sort { &numWUCont($a) <=> &numWUCont($b) } keys %{ $contig_info_ref } ) { 
    my $nt_raw = $contig_info_ref->{$id}->{'total'};
    my $nt = commify($nt_raw);

    my $frac_1 = 0;
    my $frac_2 = 0;

    if ( $contig_info_ref->{$id}->{'total'} ) { 
        $frac_1 = ( $contig_info_ref->{$id}->{'values'}->{'1'}
                   / $contig_info_ref->{$id}->{'total'} ) 
                   if ( $contig_info_ref->{$id}->{'values'}->{'1'} );
        $frac_1 = sprintf("%.02f", $frac_1);

        $frac_2 = ( $contig_info_ref->{$id}->{'values'}->{'2'}
                       / $contig_info_ref->{$id}->{'total'} )
                   if ( $contig_info_ref->{$id}->{'values'}->{'2'} );
        $frac_2 = sprintf("%.02f", $frac_2);
    }

    # Simple but efficient hack for getting mode: 
    my @modes = sort { $contig_info_ref->{$id}->{'values'}->{$b} 
                       <=> $contig_info_ref->{$id}->{'values'}->{$a} } 
                       keys %{ $contig_info_ref->{$id}->{'values'} };
    my $mode = $modes[0];

    my $sum_vals = 0;
    foreach my $i ( sort { $a <=> $b } 
                    keys %{ $contig_info_ref->{$id}->{'values'} } ) { 
        my $val = $i * $contig_info_ref->{$id}->{'values'}->{$i};
        $sum_vals += $val;
    }
    my $mean = $sum_vals / $nt_raw;
    $mean = sprintf("%.02f", $mean);

    my @output_vals = ($id, $nt, $mode, $mean, $frac_1, $frac_2);
    print_output_line(\@output_vals, $tabs, $fixed);
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

sub numWUCont { 
     my $contig = $_[0];
     if ( ( $contig !~ / Contig0*[1-9]\d* \z/xms ) 
           and   ( $contig !~ /\A Contig0 \z/xms ) ) {
         die "Malformatted input.\n";
     }
     if ($contig =~ /\A Contig0 \z/xms) {
         $contig = 0;
     }
     # Strip zeros, avoid octal-izing
     if ($contig =~ /\A Contig0*([1-9]\d*) \z/xms) { 
         $contig = $1;
     }
     return $contig;
}

sub print_output_line {
    my @_output_vals = @{ $_[0] };
    my $_tabs        = $_[1];
    my $_fixed       = $_[2];
    my $output_line  = q{};
    if (  ( $_tabs     and $_fixed     ) 
      or  ( (! $_tabs) and (! $_fixed) ) ) { 
        die "Must output either tab-delimited",
            " or fixed-width output, not both",
            " or neither!\n",
            ;
    }
    if ($_tabs) { 
        $output_line = join "\t", @_output_vals;
        print "$output_line\n";
    }
    if ($_fixed) { 
        foreach my $_val (@_output_vals) { 
            printf "%-15s", $_val;
        }
        print "\n";
    }
}

