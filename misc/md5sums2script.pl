#!/usr/bin/perl

print '#!/usr/bin/perl';
print "\n";
$report = time;

print "open REPORT, \">";
print $report;
print ".report\" ";
print "or die \"Failed to open report file\"";
print "; \n\n";

while (<>) {
    chomp($input = $_);
    if (( $input =~/^(\S+)\s+i386\/(\S+)/ ) or ( $input =~/^(\S+)\s+(\S+)/ )) { 
        print "system ";
        print '"';
        print "wget ";
        print 'ftp://ftp.toughguy.caltech.edu/pub/fedora/linux/core/updates/3/i386/';
        print $2;
        print '";';
        print "\n";

        print "system ";
        print '"';
        print "wget ";
        print 'ftp://ftp.toughguy.caltech.edu/pub/fedora/linux/core/development/i386/Fedora/RPMS/';
        print $2;
        print '";';
        print "\n";

        print "print REPORT \"Tried to get $2: \";\n";
        print "print REPORT ";
        print '`';
        print "md5sum $2 | grep $1";
        print '`;';
        print "\n";
        print "\n";
    }
}
