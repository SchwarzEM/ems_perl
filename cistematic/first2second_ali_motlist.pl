#!/usr/bin/env perl

# first2second_ali_motlist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2008.
# Purpose: Split one huge motif file into one file per contig.

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use File::Spec;

my $contig  = q{};
my $contig2fh_ref;

my $work_dir = q{};
GetOptions ( "dir=s" => \$work_dir );
$work_dir ||= q{.};
$work_dir = abs_path($work_dir);
if (! -d $work_dir ) { 
    mkdir $work_dir, 0770 
        or die "Can't create requested directory $work_dir!\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ / \A \S+ \t \d+ \t \d+ \t \S+ \t [ACGT]+ \z /xms ) { 
        warn "Could not parse input line: $input\n";
    }
    if ( $input =~ / \A (\S+) \t \d+ \t \d+ \t \S+ \t [ACGT]+ \z /xms ) { 
        $contig = $1;
        if (! exists $contig2fh_ref->{$contig} ) { 
            define_print_filehandle($contig, $contig2fh_ref);
        }
        print { $contig2fh_ref->{$contig} } "$input\n";
    }
}

close_all_contig_fhs();

sub define_print_filehandle {
    my $input_contig   = $_[0];
    my $input_c2fh_ref = $_[1];
    my $output_name = q{};
    if ($ARGV) { 
        $output_name .= $ARGV . q{.};
    }
    $output_name .= $input_contig . '.split.txt';
    $output_name = File::Spec->catfile( $work_dir, $output_name );
    $output_name = failsafe_name($output_name);
    open $input_c2fh_ref->{$input_contig}, '>', $output_name 
        or die "Can't open filehandle to $output_name: $!";
    return;
}

sub close_all_contig_fhs { 
    my $input_c2fh_ref = $_[0];
    foreach my $c ( keys %{ $input_c2fh_ref } ) { 
        close $input_c2fh_ref->{$c} 
            or die "Can't close filehandle to contig-$c-specific output: $!";
    }
    return;
}

sub failsafe_name {
    my $filename = $_[0];
    if (-e $filename) {
        my $suffix = 0;
        while (-e $filename) {
            $suffix++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix";
        }
    }
    return $filename;
}

