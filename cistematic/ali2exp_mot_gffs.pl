#!/usr/bin/env perl

# ali2uniq_mot_gffs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/23/2008.
# Purpose: from MANY Ali tables, make GFF w/ UNIQUE motif instances.
# Problem: uses long running time + lots of RAM; needs threading.

use strict;
use warnings;
use Carp;

use File::Basename;
use Scalar::Util;
use Tie::File;
use Tie::Hash;
use Time::HiRes qw ( gettimeofday tv_interval );

# Use %sites_ref to store GFFable data for site instances.
# 
# Each key is a unique site in the genome, defined 
#     as a \t-glued seq. + coord1 + coord2 triplet.
#     This scalar is a unique hash key, whose 
#     value is a *reference* to a second hash; that, 
#     in turn, is a bunch of $sourcename => 1 links.
#     So each unique site gets all its sourcenames listed,
#     making it clear if a 'unique' site was predicted >=2x.

my %sites_ref = ();
my $sitesref_on_disk = $ARGV[0] . '.SITES_REF_ON_DISK';
tie( %sites_ref, 'Tie::Hash', $sitesref_on_disk )
    or croak 'Cannot set up disk file for %sites_ref:', 
             " $!",
             ;

# Sanity check to filter out aberrantly short "motifs".
# 6-nt limit is probably arbitrary, but in real life is probably OK.
# If worthwhile, could code in command-line argument to vary this.

my $MIN_LENGTH = 6;

# All the individual motif files named in list files.
my %motfiles = ();
my $motfiles_on_disk = $ARGV[0] . '.MOTFILES_ON_DISK';
tie( %motfiles, 'Tie::Hash', $motfiles_on_disk )
    or croak 'Cannot set up disk file for %motfiles:',
             " $!",
             ;


# Track progress as the loooong program runs.
my $motlist_file_count = 0;
my $motfile_count      = 0;
my $readfile_count     = 0;
my $instance_count     = 0;

# For formatted outputs below:
my $format_string1 
    = 'Reading:  %s motif files, from %s motif list(s)';
my $format_string2 
    = 'Finished reading:  %s motif files, from %s motif list(s)';
my $format_string3  
    = 'Recorded:  %s instances from %s motif files (%.2f%% total)';  
my $format_string4 
    = 'Time elapsed:  %s seconds.  ETA: >~ %s seconds.';


# Progress reports will go to to the hard drive, 
# using Tie::File to give a small, dynamic, report file:

my @progress = ();
my $prog_rep = $ARGV[0] . '.ONGOING_PROGRESS.txt';
tie( @progress, 'Tie::File', $prog_rep ) 
    or croak "Can't manipulate progress report $prog_rep: $!";

# Read indiv. motif file names from >=1 motif-list files:
foreach my $motfile_list (@ARGV) {
    $motlist_file_count++;

    # If *any* of these fail, kill program loudly:
    open( my $MOTFILE_LIST, '<', $motfile_list )
        or croak "Can't open mot. file list $motfile_list: $!";

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
    close $MOTFILE_LIST 
        or croak "Couldn't close list of motif files $motfile_list: $!";
    push @progress, "Finished reading motif list:  $motfile_list\n";
}
    
$progress[1] = sprintf( "$format_string2\n",
                        commify($motfile_count),
                        commify($motlist_file_count),
                      );


# Use push-q{}x5 and @progress[-4,-2] to compactly update @progress.

push @progress, ( q{}, q{}, q{}, q{}, q{}, );

# Track the time elapsed during the (often long) run.
my $t0 = [gettimeofday];

foreach my $motfile (sort keys %motfiles) { 
    # Warn, don't die; no point in killing script over 1/10,000 errors.
    open(my $MOTFILE, '<', $motfile) 
        or carp "Can't open motif file $motfile: $!";

    if (Scalar::Util::openhandle $MOTFILE) { 
        $readfile_count++;
    }

    while (my $input = <$MOTFILE>) { 
        chomp $input;
        # ACGT == non-repeat; acgt == softmasked repeats, so ignore.
        # This is a strict filter; even one 'acgt' disqualifies an instance.
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
            my $gffed_seqname   = $1;
            my $first_nt        = $2;
            $first_nt++;               # Fix immediately.
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
            my $site_coords = "$gffed_seqname\t$first_nt\t$last_nt";
            # And a ref. to this site, arrayed:
            my $coord_arrayref = [$gffed_seqname, $first_nt, $last_nt];

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

            # Use a hash to have nonredundant sequences and names.
            $sites_ref{$site_coords}->{$sourcename} = 1;
            $sites_ref{$site_coords}->{'coords'} = $coord_arrayref;
            $sites_ref{$site_coords}->{'rawseq'} = $instance_rawseq;

            # Track progress so looooong runs are bearable:
            $instance_count++;
            my $eta = 0;
            if ( $readfile_count > 0 ) {
                $eta = ( ( time_elapsed($t0)
                              * ( $motfile_count - $readfile_count ) )
                            / $readfile_count );
                $eta = sprintf( "%.1f", $eta, );
            }
            $progress[-4] 
                = sprintf( "$format_string3\n", 
                           commify($instance_count), 
                           commify($readfile_count),
                           ( 100 * ( $readfile_count / $motfile_count ) ),
                         );
            $progress[-2]
                = sprintf( "$format_string4\n", 
                           commify( time_elapsed($t0) ),
                           commify($eta),
                         );
        }
    }
    # Again, don't let 1/10,000 failures kill whole script:
    close $MOTFILE 
        or carp "Couldn't close individual motif file $motfile: $!";
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

# Sort coordinates by ascent through genome:
foreach my $site ( grep { $_ =~ /\S/xms } 
                   sort { $sites_ref{$a}->{'coords'}->[0] 
                              cmp $sites_ref{$b}->{'coords'}->[0] }
                   sort { $sites_ref{$a}->{'coords'}->[1] 
                              <=> $sites_ref{$b}->{'coords'}->[1] }
                   sort { $sites_ref{$a}->{'coords'}->[2] 
                              <=> $sites_ref{$b}->{'coords'}->[2] }
                   keys %sites_ref ) { 

    # Get list of any multiple hits on a unique sequence.
    my @motif_names = grep { $_ ne 'coords'} 
                      grep { $_ ne 'rawseq'}
                      sort keys %{ $sites_ref{$site} } ;
    my $motif_namelist = join q{|}, @motif_names;
    my $rawseq = $sites_ref{$site}->{'rawseq'};

    print "$site",           # '\w+\t\d+\t\d+' coord. string.
          "\t",
          "$motif_namelist",
          "\t",
          $rawseq,
          "\n",
    ;

}


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $text = reverse $_[0];
    $text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $text;
}

sub time_elapsed {
    my $t0 = $_[0];
    my $time_elapsed 
        = sprintf ( "%.1f", tv_interval($t0, [gettimeofday]), );
    return $time_elapsed;
}

