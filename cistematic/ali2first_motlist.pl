#!/usr/bin/env perl

# ali2uniq_mot_gffs_slow.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/31/2008.
# Purpose: given list of MANY Ali tables, make 1 flatfile per contig/chromosome, sorted/nonredund. wrt each motif file.
# Have sorted nonredundancy of lines read from a motif file, but otherwise do not
#     try to sort or consolidate data.

use strict;
use warnings;

use Carp;
use File::Basename;
use File::Path;
use File::Spec;
use Tie::File;
use Time::HiRes qw ( gettimeofday tv_interval );

# Have unambiguous, permanently usable stem for all file names:

my $job_name = $ARGV[0];
if ( (! $job_name ) or ( $job_name eq '-' ) or ( $job_name !~ /\S/xms ) ) { 
    $job_name = 'stdin';
}

### TRACK PROGRESS OF RUNS (which can be long): ###

my $date  = join('.', &get_local_date());
$job_name = $job_name . q{.} . $date;

# Name array and disk file for dynamic progress reporting:
my @progress = ();
my $prog_rep      = $job_name . '.PROGRESS.txt';
$prog_rep         = failsafe_name($prog_rep);
tie( @progress, 'Tie::File', $prog_rep )
    or croak "Can't manipulate progress report $prog_rep: $!";

# For formatted progress reports below.
my $motlist_file_count = 0;
my $motfile_count      = 0;
my $readfile_count     = 0;
my $instance_count     = 0;

### STORE LIST of INDIV. MOTIF FILES (from >=1 motif-list files). ###

# Store all named individual motif files here:
my %motfiles = ();

foreach my $motfile_list (@ARGV) {
    open my $MOTFILE_LIST, '<', $motfile_list 
        or croak "Can't open motif-file list $motfile_list: $!";

    $motlist_file_count++;

    while ( my $input = <$MOTFILE_LIST> ) { 
        chomp $input;
        $motfiles{$input} = 1;
        $motfile_count++;
        my $format1 = 'Reading:  %s motif files, from %s motif list(s)';
        $progress[1] 
            = sprintf( "$format1\n",
                       commify($motfile_count),
                       commify($motlist_file_count),
                     );
    }

    close $MOTFILE_LIST 
        or croak "Couldn't close filehandle from $motfile_list: $!";

    push @progress, "Finished reading motif list:  $motfile_list\n";
}


### EXPAND ONGOING PROGRESS REPORT: ###

my $format2 = 'Finished reading:  %s motif files, from %s motif list(s)';
$progress[1] = sprintf( "$format2\n",
                        commify($motfile_count),
                        commify($motlist_file_count),
                      );


# Use push-q{}x5 and @progress[-4,-2] to compactly update @progress.
push @progress, ( q{}, q{}, q{}, q{}, q{}, );

# Track the time elapsed during the (often long) run.
my $t0 = [gettimeofday];

# Use $MIN_LENGTH to filter out absurdly short "motifs".
# A 6-nt limit is probably arbitrary, but in real life is probably OK.  
# If worthwhile, one could code in a command-line argument to vary this.

my $MIN_LENGTH = 6;

### EXTRACT DATA FROM (perhaps very many) INDIV. MOTIF FILES: ###

foreach my $motfile (sort keys %motfiles) { 
    # Build up and print nonredundant, sorted data for each motfile.
    my %motdata = ();

    # Warn, don't die; no point in killing the script over 1/10,000 errors!
    open my $MOTFILE, '<', $motfile 
        or carp "Can't open motif file $motfile: $!";
    $readfile_count++;

    while (my $input = <$MOTFILE>) { 
        chomp $input;
        # ACGT == non-repeat; acgt == softmasked repeats, so ignore.
        # This is a strict filter:
        #   just one lowercase 'a', 'c', 'g', or 't' disqualifies an instance.
        if ( $input =~ / (?:[^\t]*\t){3} 

                         # $1 -> name of seq.
                         ( \S+ ) 
                         \t

                         # $2 -> start nt, with one-off-low error:
                         ( \d+ ) 
                         \t

                         # $3 -> site orientation, either 'F' or 'R'
                         ( F | R ) 
                         \t
                         (?:[^\t]*\t){6}

                         # $4 -> raw sequence (to get site length):
                         ( [ACGT]{$MIN_LENGTH,} ) 
                         \t 

                       /xms ) { 

            # Extract out GFF-able data from an Ali table.
            my $gff_seqname     = $1;
            my $first_nt        = $2;
            $first_nt++;               # Fixed immediately.
            my $orientation     = $3;
            my $instance_rawseq = $4;

            # Work out GFF-able coordinates:
            my $site_reach   = ( length($instance_rawseq) - 1);

            # From $first_nt, far end of site is $site_reach nt away.
            my $last_nt      = q{};
            if ( $orientation eq 'F' ) {
                 $last_nt = $first_nt + $site_reach;
            } 
            if ( $orientation eq 'R' ) { 
                 $last_nt = $first_nt;
                 $first_nt = $first_nt - $site_reach;
            }

            # Get name of motif from its filename.
            my $sourcename = basename($motfile); 

            # Trim off '-ce.\d+.tab' suffix:
            $sourcename =~ s/\-ce\.\d+\.tab\z//xms;

            # Use generic name for unnamed motifs or STDIN.
            if (  ( $sourcename eq q{-}  ) 
               or ( $sourcename !~ /\w/xms ) ) {
                # Kludge: need elegant sorting of this to end of list.
                $sourcename = 'unspec_motif';
            }

            # Define a single data string with '\t' splits between fields:            
            my $data_line = $gff_seqname
                            . "\t"
                            . $first_nt
                            . "\t"
                            . $last_nt
                            . "\t"
                            . $sourcename
                            . "\t"
                            . $instance_rawseq
                            ;
            $motdata{$data_line} = 1;

            # Track progress so looooong runs are bearable:
            $instance_count++;
            my $eta = 0;
            if ( $readfile_count > 0 ) {
                $eta = ( ( time_elapsed($t0)
                              * ( $motfile_count - $readfile_count ) )
                            / $readfile_count );
                $eta = sprintf( "%.1f", $eta, );
            }
            my $format3 = 'Recorded:  %s instances from %s motif files (%.2f%% total)';
            $progress[-4] 
                = sprintf( "$format3\n", 
                           commify($instance_count), 
                           commify($readfile_count),
                           ( 100 * ( $readfile_count / $motfile_count ) ),
                         );
            my $format4 = 'Time elapsed:  %s seconds.  ETA: >~ %s seconds.';
            $progress[-2]
                = sprintf( "$format4\n", 
                           commify( time_elapsed($t0) ),
                           commify($eta),
                         );
        }
    }
    # Again, don't let 1/10,000 failures kill whole script:
    close $MOTFILE 
        or carp "Couldn't close individual motif file $motfile: $!";
    
    # Then emit print of nonredundant sorted contents of the file:

    foreach my $output ( sort keys %motdata ) { 
        print "$output\n";
    }
}

### UPDATE ONGOING PROGRESS REPORT: ###

$progress[-4] 
    = sprintf( "Extracted: %s instances out of %s motif files\n", 
               commify($instance_count), 
               commify($motfile_count),
             );

$progress[-2]
    = sprintf( "Time elapsed: %s seconds.\n",
               commify( time_elapsed($t0) ),
             );

### SUBROUTINES: ###

sub commify { 
    my $text = reverse $_[0];
    $text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $text;
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

sub time_elapsed {
    my $t0 = $_[0];
    my $time_elapsed 
        = sprintf ( "%.1f", tv_interval($t0, [gettimeofday]), );
    return $time_elapsed;
}

