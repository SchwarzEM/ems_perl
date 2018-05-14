#!/usr/bin/perl

# when clustering NR, it core-dumped after writting the database
# before output the .clstr file
# This script is to generate this .clstr file from .bak.clstr

my %cluster = ();
while($ll=<>) {
  chop($ll);
  my ($i, $v) = split(/\t/,$ll);
  if ($i < 0 ) {$i = -$i-1;}

  if (not defined($cluster{$i})) {
    $cluster{$i} = [$v];
  }
  else {
    push(@{$cluster{$i}},$v);
  }
}

my @c = sort {$a <=> $b} keys %cluster;

for ($i=0; $i<=$#c; $i++) {
  print ">Cluster $i\n";

  $v = $c[$i];
  for ($j=0; $j<=$#{$cluster{$v}}; $j++) {
    print "$j\t$cluster{$v}->[$j]\n";
  }
}
