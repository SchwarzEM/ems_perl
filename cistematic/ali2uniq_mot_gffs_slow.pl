#!/usr/bin/env perl

# ali2uniq_mot_gffs_slow.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/31/2008.
# Purpose: from list of MANY Ali tables, make GFF w/ UNIQUE motif instances.
# Problems: EXTREMELY long runtime; not really usable for big-genome work.

use strict;
use warnings;

use Carp;
use DBM::Deep;
use File::Basename;
use File::Path;
use File::Spec;
use Tie::File;
use Time::HiRes qw ( gettimeofday tv_interval );

# Have unambiguous, permanently usable stem for all file names:

my $job_name = $ARGV[0];
if ( ( $job_name eq '-' ) or ( $job_name !~ /\S/xms ) ) { 
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

# Name directory where BIG hashref files (1 per contig or chromosome) will go.
my $job_dir = $job_name . '.dir';
$job_dir    = failsafe_name($job_dir);
mkdir $job_dir, 0750;

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

### INITIALIZE VARIABLES FOR STORING GFF-able SITE INSTANCES: ###

# $sites_ref is a storage hashREF, NOT itself tied to DBM::Deep.
#     Its keys are the sequences listed nonredundantly in @sequences.
#     And its *values* are DBM::Deep objects! i.e., blessed hashrefs!:
#         $sites_ref->{$seqN} = DBM::Deep->new( $workdir/($seqN .'.db') );
#     This wackery should let me do arbitrarily many DBMs ad lib.
# 
# Each $sites_ref->{$seqN} value is itself a reference to a hash, 
#     for which each key is a unique site in the genome, defined as
#     a \t-glued seq. + coord1 + coord2 triplet.
# 
# This key refers, in turn, to a hashref; that hash
#     is a bunch of $sourcename => 1 links.
# 
# So:
#     each sequence gets a DBM file
#     each DBM file stores a hash
#     whose keys are unique genome sites
#     linked to nonredundantly listed sourcenames
#     making it clear if a 'unique' site was predicted >=2x.

my $sites_ref;

# Use $MIN_LENGTH to filter out absurdly short "motifs".
# A 6-nt limit is probably arbitrary, but in real life is probably OK.  
# If worthwhile, one could code in a command-line argument to vary this.

my $MIN_LENGTH = 6;

### EXTRACT DATA FROM (perhaps very many) INDIV. MOTIF FILES: ###

foreach my $motfile (sort keys %motfiles) { 
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

            # Define a triplet site coordinate string:
            my $site_coords = "$gff_seqname\t$first_nt\t$last_nt";

            # And a ref. to this site, arrayed:
            my $coord_arrayref = [$gff_seqname, $first_nt, $last_nt];

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

            # Use a seq-spec. hash:
            #    partly to have nonredundant sequences and names.
            #    partly to break up VAST predictions into <2-GB files.
            # 
            # N.B.: this relies on $sites_ref->{"$job_dir/$gff_seqname"} 
            #     having actually been made a DBM::Deep-blessed hashref.

            ensure_DBM_Deep_hashref($sites_ref, $gff_seqname, $job_dir);

            $sites_ref->{$gff_seqname}->{$site_coords}->{$sourcename} 
                = 1;
            $sites_ref->{$gff_seqname}->{$site_coords}->{'coords'} 
                = $coord_arrayref;
            $sites_ref->{$gff_seqname}->{$site_coords}->{'rawseq'} 
                = $instance_rawseq;

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


### SORT COORDINATES BY ASCENT THROUGH GENOME, then PRINT THEM ALL: ###

# This is apparently where things failed in a big run.
# Should probably code in a Schwartzian transform, if I can.

my @sequences = grep { $_ =~ /\S/xms } sort keys %{ $sites_ref };

foreach my $seq (@sequences) { 
    my @sites_list = ();
    my $sites_array = $job_name . '.sites_array';
    $sites_array = File::Spec->catfile( $job_dir, $sites_array );
    $sites_array = failsafe_name($sites_array);

    tie( @sites_list, 'Tie::File', $sites_array )
        or croak "Can't load sites array into $sites_array: $!";

    @sites_list = grep { $_ =~ /\S/xms }
                  sort { $sites_ref->{$seq}->{$a}->{'coords'}->[0]
                             cmp $sites_ref->{$seq}->{$b}->{'coords'}->[0] }
                  sort { $sites_ref->{$seq}->{$a}->{'coords'}->[1]
                             <=> $sites_ref->{$seq}->{$b}->{'coords'}->[1] }
                  sort { $sites_ref->{$seq}->{$a}->{'coords'}->[2]
                             <=> $sites_ref->{$seq}->{$b}->{'coords'}->[2] }
                  keys %{ $sites_ref->{$seq} };

    foreach my $site (@sites_list) { 
        # Get list of any multiple hits on a unique sequence.
        my @motif_names = grep { $_ ne 'coords'} 
                          grep { $_ ne 'rawseq'}
                          sort keys %{ $sites_ref->{$seq}->{$site} } ;
        my $motif_namelist = join q{|}, @motif_names;
        my $rawseq = $sites_ref->{$seq}->{$site}->{'rawseq'};

        print "$site",           # '\w+\t\d+\t\d+' coord. string.
              "\t",
              "$motif_namelist",
              "\t",
              $rawseq,
              "\n",
              ;
    }
}


### CLEAN UP and QUIT: ###

$sites_ref = q{};

# Left commented out for now, until I've had one successful big run:
#
# rmtree($job_dir, 0, 1);


### SUBROUTINES: ###

sub commify { 
    my $text = reverse $_[0];
    $text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $text;
}

sub ensure_DBM_Deep_hashref {
    my ($top_hashref, $seqname, $directory) = @_;

    # Defend against lousy seqnames:
    if ( $seqname =~ / \A \. /xms ) { 
        warn "'$seqname' has a leading '.',",
             " and thus makes files invisible when",
             " used as a prefix!\n",
             ;
    }
    if ( $seqname !~ /\S/xms ) { 
        die "'$seqname' is not a usable sequence name.\n";
    }

    # If the hash entry exists, do nothing quietly:
    if ( exists $top_hashref->{$seqname} ) { 
        return;
    }

    # If the hash entry does not exist, make a DBM::Deep object.
    if (! exists $top_hashref->{$seqname} ) { 
        my $tiefile_path = $seqname . q{.} . 'db';
        # Use a working directory if it's defined and exists.
        if ( ( $directory ) and ( -e $directory ) ) { 
            $tiefile_path 
                = File::Spec->catfile( $directory, $tiefile_path );
        }
        $top_hashref->{$seqname} = DBM::Deep->new( $tiefile_path ) 
            or die "$tiefile_path unsuccessfully DBM::Deep-blessed: $!";
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

