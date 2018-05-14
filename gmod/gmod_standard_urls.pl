#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use Bio::GMOD::StandardURLs;

$|++;

my $requested_mod  = shift;

$requested_mod or die "Usage: gmod_standard_urls.pl [mod] eg WormBase, FlyBase, SGD)";

my $mod = Bio::GMOD::StandardURLs->new(-mod=>$requested_mod);
my $available_species = $mod->available_species(-expanded=>1);

print "$requested_mod provides data for the following species using GMOD common URLs\n";
print map { "\t$_\n" } keys %$available_species;
print "\n";

foreach my $binomial (keys %$available_species) {
  print "$binomial\n";
  print '-' x length ($binomial) . "\n";
   my @releases = $mod->releases(-species=>$binomial,-expanded=>1);
   foreach (@releases) {
       my ($release,$date,$availability) = @$_;
       print join("\t",$release,$date,$availability),"\n";
       next if ($availability eq 'no');
       my $datasets = $mod->datasets(-species=>$binomial,-release=>$release);
       foreach (keys %$datasets) {
            print "\t$_\t$datasets->{$_}\n";
       }
   }
  print "\n";
}


__END__

=pod

=head1 NAME

gmod_fetch_features.pl - 
gmod_update_installation.pl - Maintain a MOD installation

=head1 USAGE

This script provides a convenient mechanism to maintain a MOD
installation.  It should be excecuted with super user privileges.

  $ gmod_update_installation.pl [options]

=head1 OPTIONS

The following options are generically available for any MOD (default
values in parenthesis):

 MOD:
 --mod       One of WormBase, FlyBase, SGD, etc

 Versions:
 --sync_to   [live || dev] Sync to the current live or development version (live)
 --force     [boolean] Force an update to the live or development version as appropriate (false)
 --version   Update to the provided version (the current live version)

 System paths:
 --tmp       Full path to the temporary directory to hold downloads (/usr/local/gmod/tmp)

 Miscellaneous:
 --purge     [boolean] Purge the tmp download folder following upgrade (false)
 --help      Display this message

Due to the wide variety of installation paths and MOD structures, each
MOD may offer specialized options.  These can be provided as
"--option_name OPTION" which will be passed directly to the
Bio::GMOD::Update::"MOD" object's update() method.  For example, a
typical command to maintain a WormBase installation looks like:

 % gmod_update_installation.pl --analyze_logs --mysql_path /usr/local/mysql/data

For a full description of all available system paths and update
options for your particular MOD, see L<Bio::GMOD::Adaptor> and
L<GMOD::Adaptor::your_mod>.

=head1 Running under cron

You may wish to run this script under cron to ensure that your
installation is always up-to-date.  For my personal installation of
WormBase, I use the following settings:

0 2 * * * /usr/local/bin/gmod_update_installation.pl --sync_to dev

This will check for and install a new version if present at 2 AM in
the morning.

I keep my installation in sync with the development version.  You will
want to use the more stable live version, which you can specify using
"--sync_to live" or by simply leaving off the "--sync_to" option
altogether.

A suggested crontab entry for a simple local installation is:

  gmod_update_intallation.pl --sync_to live --purge 1

A suggested crontab entry for official WormBase mirror sites is:

  gmod_update_intallation.pl --sync_to live --purge 1 --analyze_logs 1

=head1 SEE ALSO

L<Bio::GMOD>, L<Bio::GMOD::Update>

=head1 AUTHOR

Todd Harris <harris@cshl.edu>.

Copyright (c) 2003-2005 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut


