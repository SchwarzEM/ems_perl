#!/usr/bin/env perl

use strict;
use warnings ;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

my @infiles;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: total_pme_expected_count_15feb2012.pl\n",
        "    --infile|-i   <input stream/files>\n",
        "    --help|-h     [print this message]\n",
        ;
}

foreach my $infile (@infiles) {
    my $basename = basename $infile;
    $basename =~ s/\.genes\.results//;
    my $total_pme_expected_count = 0;
    open my $INFILE, '<', $infile or die "Can't open input file $infile\n";
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input =~ /\A ([^\t]* \t){7} (\S+) \t/xms ) { 
            my $pme_expected_count = $1;
            if ( looks_like_number($pme_expected_count) ) { 
                $total_pme_expected_count += $pme_expected_count;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile\n";
    my $int_total_pme_expected_count = int($total_pme_expected_count);
    $int_total_pme_expected_count = commify($int_total_pme_expected_count);
    print "$basename\t$int_total_pme_expected_count\n";
}

sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

