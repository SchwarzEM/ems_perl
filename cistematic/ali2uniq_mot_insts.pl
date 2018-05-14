#!/usr/bin/env perl

# ali2uniq_mot_insts.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/27/2008.
# Purpose: from MANY Ali tables, make FASTA w/ UNIQUE motif instances.
# Problem: uses long running time + lots of RAM; needs threading.

use strict;
use warnings;
use File::Basename;
use Scalar::Util;
use Tie::File;
use Time::HiRes qw ( gettimeofday tv_interval );

# Use to store unique seqs. and their origins.
# 
# Each key is a unique nucleotide sequence, whose 
#     value is a *reference* to a hash; the hash, 
#     in turn, is a bunch of $sourcename => 1 links.
#     So each seq. gets all its sourcenames listed.
# 
my %seqs2names_ref = ();

# Sanity check to filter out aberrantly short "motifs".
# If worthwhile, could code in command-line argument to vary this.
my $MIN_LENGTH = 6;

my %motfiles = ();

# Track progress as the long program runs.
my $motlist_file_count = 0;
my $motfile_count      = 0;
my $readfile_count     = 0;
my $instance_count     = 0;

my $format_string1 
    = 'Reading:  %s motif files, from %s motif list(s)';
my $format_string2 
    = 'Finished reading:  %s motif files, from %s motif list(s)';
my $format_string3  
    = 'Recorded:  %s instances from %s motif files (%.2f%% total)';  
my $format_string4 
    = 'Time elapsed:  %s seconds.  ETA: >~ %s seconds.';

# Have unambiguous part of prog-file names, to avoid overwriting by diff. runs.
my $date = join('.', &get_local_date());

# To avoid annoyingly long outputs, remember short "progress" string ... 
my @progress = ();
my $prog_rep = $ARGV[0] . q{.} . $date . '.PROGRESS.txt';

# ... and use Tie::File to allow transparent revisions w/o massive file:
tie(@progress, "Tie::File", $prog_rep) 
    or die "Can't manipulate progress report $prog_rep: $!";

foreach my $motfile_list (@ARGV) {
    $motlist_file_count++;
    open(my $MOTFILE_LIST, "<", $motfile_list) 
        or die "Can't open mot. file list $motfile_list: $!";
    while (my $input = <$MOTFILE_LIST>) { 
        chomp $input;
        $motfiles{$input} = 1;
        $motfile_count++;
        $progress[1] 
            = sprintf( "$format_string1\n",
                       commify($motfile_count),
                       commify($motlist_file_count),
                     );
    }
    close $MOTFILE_LIST;
    push @progress, "Finished reading motif list:  $motfile_list\n";
}
    
$progress[1] = sprintf( "$format_string2\n",
                        commify($motfile_count),
                        commify($motlist_file_count),
                      );


# Use push-q{}x5 and @progress[-4,-2] to compactly update @progress.

push @progress, ( q{}, q{}, q{}, q{}, q{});

# Track the time elapsed during the (often long) run.
my $t0 = [gettimeofday];

foreach my $motfile (sort keys %motfiles) { 
    # Warn, don't die; no point in killing script over 1/10,000 errors.
    open(my $MOTFILE, "<", $motfile) 
        or warn "Can't open motif file $motfile: $!";

    if (Scalar::Util::openhandle $MOTFILE) { $
        readfile_count++;
    }

    while (my $input = <$MOTFILE>) { 
        chomp $input;
        # ACGT == non-repeat; acgt == softmasked repeats, so ignore.
        if ( $input =~ / (?:[^\t]*\t){12} 
                         ( [ACGT]{$MIN_LENGTH,} ) 
                         \t 
                       /xms ) { 

            # Extract out raw sequence from an Ali table.
            my $sequence = $1;

            # Get name of motif from its filename.
            my $sourcename = basename($motfile); 
            # Trim off '-ce.\d+.tab' suffix:
            $sourcename =~ s/\-ce\.\d+\.tab\z//xms;

            # Use generic name for unnamed motifs or STDIN.
            if (  ( $sourcename eq "-"  ) 
               or ( $sourcename !~ /\w/ ) ) {
                # Kludge: need elegant sorting of this to end of list.
                $sourcename = "unspec_motif";
            }

            # Use a hash to have nonredundant sequences and names.
            $seqs2names_ref{$sequence}->{$sourcename} = 1;
            $instance_count++;
            $progress[-4] 
                = sprintf( "$format_string3\n", 
                           commify($instance_count), 
                           commify($readfile_count),
                           ( 100 * ( $readfile_count / $motfile_count ) ),
                         );
            my $eta = 0;
            if ( $readfile_count > 0 ) { 
                $eta = ( ( time_elapsed($t0) 
                              * ( $motfile_count - $readfile_count ) ) 
                            / $readfile_count );
                $eta = sprintf( "%.1f", $eta, );
            }
            $progress[-2]
                = sprintf( "$format_string4\n", 
                           commify( time_elapsed($t0) ),
                           commify($eta),
                         );
        }
    }
    close $MOTFILE;
}

# Record next increment:
$progress[-4] 
    = sprintf( "Extracted: %s instances out of %s motif files\n", 
               commify($instance_count), 
               commify($motfile_count),
             );

$progress[-2]
    = sprintf( "Time elapsed: %s seconds.\n",
               commify( time_elapsed($t0) ),
             );

my $i = 1;
foreach my $seq (grep { $_ =~ /\S/ } sort keys %seqs2names_ref ) { 

    # Get list of any multiple hits on a unique sequence.
    my @seqnames = sort keys %{ $seqs2names_ref{$seq} } ;
    my $seqnamelist = join "|", @seqnames;

    # All motif nos. hard-coded nine leading 0s; should code dynamically.
    my $j = sprintf("%09d", $i);

    print ">",
          "$seqnames[0]",
          "_",
          "$j\t",
          "$seqnamelist\n";
    ;
    print "$seq\n";

    $i++;
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

sub time_elapsed {
    my $_t0 = $_[0];
    my $_time_elapsed 
        = sprintf ( "%.1f", tv_interval($_t0, [gettimeofday]), );
    return $_time_elapsed;
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

