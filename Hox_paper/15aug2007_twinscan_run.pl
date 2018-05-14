#!/usr/bin/perl

# 15aug2007_twinscan_run.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/15/2007.
# Purpose: yet another quick-and-dirty gene prediction set for PB2801 supercontigs.

use strict;
use warnings;
use File::Spec;

# Example of how to run:
# ./15aug2007_twinscan_run.pl                             \
#    /usr/local/src/Twinscan_3.5_src                      \
#    bin/iscan                                            \
#    parameters/worm_iscan_est-5584-genes-6-28-2005.zhmm  \
#    /home/schwarz/taygeta/science/CB5161/2006.10.CB5161.PS1010seqs/CB5161    # or .../PS1010

my $date = join('.', &get_local_date());

my $twinscan_dir = $ARGV[0];
my $tscan_bin    = $ARGV[1];
my $tscan_params = $ARGV[2];
my $seq_dir      = $ARGV[3];

$tscan_bin = File::Spec->catfile($twinscan_dir, $tscan_bin);
if (! ((-e $tscan_bin) and (-x $tscan_bin)) ) {
    die "Can't execute $tscan_bin: $!";
}
$tscan_params = File::Spec->catfile($twinscan_dir, $tscan_params);
if (! ((-e $tscan_params) and (-r $tscan_params)) ) {
    die "Can't read $tscan_params: $!";
}
if (! ((-e $seq_dir) and (-d $seq_dir )) ) {
    die "Can't use $seq_dir as directory: $!";
}
opendir my $SEQ_DIR, $seq_dir or die "Can't open directory $seq_dir: $!";
if (-e "$date.gff") {
    system "rm $date.gff";
}
my @seq_files = sort ( grep { $_ !~ /^\./ } (readdir $SEQ_DIR) );

foreach my $seq_file (@seq_files) { 
    $seq_file = File::Spec->catfile($seq_dir, $seq_file);
    system "$tscan_bin $tscan_params $seq_file 1>>$date.gff 2>>errors.$date.txt";
}

sub get_local_date { 
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

