#!/usr/bin/env perl

# parse_nyc_covid19_29mar2020.pl -- Erich Schwarz <ems394@cornell.edu>, 3/29/2020.
# Purpose: automatic download, filtering, human-readable formatting, and printout of Covid-19 numbers for New York City.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $date   = q{};
my $cases  = 0;
my $deaths = 0;

my $prev_cases  = 0;
my $prev_deaths = 0;

my $new_cases  = 0;
my $new_deaths = 0;

my $daily_new_case_ratio  = 0;
my $daily_new_death_ratio = 0;

my $death_to_case_ratio   = 0;

my $infile = q{};

my $header =   "   Date         "
             . "Cases    "
             . "New_cases   "
             . "Case_ratio   "
             . "Deaths   "
             . "New_deaths   "
             . "Death_ratio   "
             . "Death/Cases"
             . "\n"
             ;

my $help;

GetOptions ( 'infile=s' => \$infile,
             'help'     => \$help,   );

# If an existing infile has not been named), download a new CSV file and timestamp its name down to one second.
if (! $infile ) {
    my $time = join( '.', &get_local_date() );
    $infile  = "$time.us-counties.csv";
    system "curl --silent --output $infile https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv ;";
}

if ( $help or (! -e $infile) ) { 
    die "Format: parse_nyc_covid19_29mar2020.pl\n",
        "    --infile|-i     [pre-existing input CSV file to parse; overrides default of downloading CSV file from NYT]]\n",
        "    --help|-h       <print this message>\n",
        "    [prints to STDOUT; redirect to a file with '>']\n";
        ;
}

my %date_comment = (
    '2020-03-01' => '   # Sun. 3/1, start of available records:',
    '2020-03-14' => '   # Sat. 3/14, first recorded Covid-19 death:',
    '2020-03-16' => '   # Mon. 3/16, De Blasio defends his going to the gym:',
    '2020-03-22' => '   # Sun. 3/22, NYC goes into lockdown ordered by Cuomo:',
    '2020-03-27' => '   # Fri. 3/27, trials of hydroxychloroquine and azithromycin begin on ~1,100 patients:',
);

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;

    # Strip out noise from each line of the alleged CSV data file:
    $input =~ s/\A\s+//;
    $input =~ s/\<.+?\>//g;

    if (     ( $input =~ /\A [^,]+ [,] New [ ] York [ ] City [,] /xms ) 
         and ( $input =~ /\A ([^,]*) [,] 
                              [^,]*  [,] 
                              [^,]*  [,] 
                              [^,]*  [,] 
                             ([^,]*) [,] 
                             ([^,]*)\z/xms ) 
       ) { 
        $date   = $1;
        $cases  = $2;
        $deaths = $3;

        $new_cases  = $cases - $prev_cases ;
        $new_deaths = $deaths - $prev_deaths ;
        if ($prev_cases) {
            $daily_new_case_ratio = ($cases/$prev_cases);
            $daily_new_case_ratio = sprintf("%.2f", $daily_new_case_ratio);
        }
        if ($prev_deaths) {
             $daily_new_death_ratio = ($deaths/$prev_deaths);
             $daily_new_death_ratio = sprintf("%.2f", $daily_new_death_ratio);
        }

        if ($cases) {
             $death_to_case_ratio = ($deaths/$cases);
             $death_to_case_ratio = (100 * $death_to_case_ratio);
             $death_to_case_ratio = sprintf("%.2f", $death_to_case_ratio);
             $death_to_case_ratio = $death_to_case_ratio . '%';
        }


        print "\n$header\n" if $header;
        $header = q{};

        my $print_cases     = &commify( $cases );
        my $print_new_cases = &commify( $new_cases );

        my $print_deaths     = &commify( $deaths );
        my $print_new_deaths = &commify( $new_deaths );
 
        if ( exists $date_comment{$date} ) {
            print "\n" unless $date eq '2020-03-01';
            print "$date_comment{$date}\n";
        }

        $date                  = sprintf("%-13s",$date);
        $print_cases           = sprintf("%-9s",$print_cases);
        $print_new_cases       = sprintf("%-12s",$print_new_cases);
        $daily_new_case_ratio  = sprintf("%-13s",$daily_new_case_ratio);
        $print_deaths          = sprintf("%-9s",$print_deaths);
        $print_new_deaths      = sprintf("%-13s",$print_new_deaths);
        $daily_new_death_ratio = sprintf("%-14s",$daily_new_death_ratio);
        $death_to_case_ratio   = sprintf("%-14s",$death_to_case_ratio);

        print "   $date",
              "$print_cases",
              "$print_new_cases",
              "$daily_new_case_ratio",
              "$print_deaths",
              "$print_new_deaths",
              "$daily_new_death_ratio",
              "$death_to_case_ratio",
              "\n",
              ;

        $prev_cases  = $cases;
        $prev_deaths = $deaths ;
        $new_cases   = 0;
        $new_deaths  = 0;

        $daily_new_case_ratio  = 0;
        $daily_new_death_ratio = 0;

        $print_cases      = 0;
        $print_new_cases  = 0;
        $print_deaths     = 0;
        $print_new_deaths = 0;
    }
    # Enforce correct parsing of any line containing "New York City" in it:
    if (     ( $input =~ /\A New [ ] York [ ] City /xms )  
         and ( $input !~ /\A ([^,]*) [,]
                              [^,]*  [,]
                              [^,]*  [,]
                              [^,]*  [,]
                             ([^,]*) [,]
                             ([^,]*)\z/xms ) 
       ) {
        die "Cannot parse: $input\n";
    }
}
print "\n";
close $INFILE;

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

sub get_local_date {
    my @_ltime = localtime;
    my @_ldate = ( (sprintf ("%04u", ($_ltime[5] + 1900)) ),     # $year
                   (sprintf ("%02u", ($_ltime[4] + 1))    ),     # $mon
                   (sprintf ("%02u", ($_ltime[3] + 0))    ),     # $mday
                   (sprintf ("%02u", ($_ltime[2] + 0))    ),     # $hour
                   (sprintf ("%02u", ($_ltime[1] + 0))    ),     # $min
                   (sprintf ("%02u", ($_ltime[0] + 0))    ), );  # $sec
    return @_ldate;
}

