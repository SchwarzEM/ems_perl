#!/usr/bin/env perl

use strict;
use warnings;

while (my $main_script = <>) { 
    chomp $main_script;
    my $backup_script = $main_script . '.backup';
    my $inline_perl_text_patch = q{ perl -ne ' s/24:00:00/08:00:00/; print ' };
    system "    mv -i $main_script $backup_script ;\n";
    system "    cat $backup_script | $inline_perl_text_patch > $main_script ;\n";
}

