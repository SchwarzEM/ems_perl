#!/usr/bin/env perl

use strict;
use warnings;

my %lib2class = ( '10862' => '300 nt',
                  '10946' => '300 nt',
                  '10kb'  => '10 kb',
                  '11318' => '300 nt',
                  '11800' => '500 nt',
                  '2kb'   => '2 kb',
                  '5kb'   => '5 kb', );

my %read_prefix = ( 'pe' => '2x',
                    'se' => q{}, );

my @infiles     = @ARGV;
my $genome_size = 315_000_000;
my $total_nt    = 0;
my $data_ref;

my $header = "\nInsert size\tLibrary ID\tEnd type\tReads\tNt\tCoverage\n";

foreach my $infile (@infiles) { 
    if ( $infile =~ / \A done\.ok\/Hco_v4_genDNA\.(\S+)\.SOAPdn\.NONcontam_reads\.((?:p|s)e)\.fa\.nt\.count\.txt \z/xms ) { 
        my $lib_name = $1;
        my $end_type = $2;
        my $nt       = 0;
        my $reads    = 0;

        open my $INFILE, '<', $infile or die "Can't open input file: $infile\n";

        while (my $input = <$INFILE>) { 
            chomp $input;
            if ( $input =~ /\A Total \s+ nt: \s+ (\S+) \s* \z/xms ) {
                $nt = $1;
            }
            elsif ( $input =~ /\A Scaffolds: \s+ (\S+) \s* \z/xms ) {
                $reads = $1;
            }
        }

        $nt       =~ s/,//g;
        $reads    =~ s/,//g;
        $total_nt += $nt;

        if (! exists $lib2class{$lib_name}) { 
            die "Can't parse library name $lib_name\n";
        }

        my $lib_class = $lib2class{$lib_name};
        my $coverage  = ($nt / $genome_size);
        $coverage = sprintf "%.1f", $coverage;

        print $header if $header;
        $header = q{};

        if ( $end_type eq 'pe' ) {
            $reads = ($reads / 2);
        }

        $reads = commify($reads);
        $nt    = commify($nt);

        print "$lib_class\t$lib_name\t$end_type\t$read_prefix{$end_type}$reads\t$nt\t$coverage\n";

        close $INFILE or die "Can't close filehandle to input file: $infile\n";
    }
    else { 
        die "Can't parse this particular file: $infile\n";
    }
}

my $total_coverage = ($total_nt/ $genome_size);
$total_coverage    = sprintf "%.1f", $total_coverage;

$total_nt = commify($total_nt);

print "\n";
print "Total: $total_nt nt; $total_coverage", "x coverage\n", ;
print "\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

