#! /usr/bin/perl -w

use strict;

system "cp dist/sTASSEL.jar ../tassel4.0_standalone/sTASSEL.jar";


my $versionComment4 = `./run_pipeline_non_verbose.pl -versionComment`;
$versionComment4 =~ s/\s+$//;
my $versionTag4 = `./run_pipeline_non_verbose.pl -versionTag`;
$versionTag4 =~ s/\s+$//;

if (!defined $versionComment4 || length $versionComment4 == 0) {
   print "versionComment4 is not defined\n";
   exit;
}

if (!defined $versionTag4 || length $versionTag4 == 0) {
   print "versionTag4 is not defined\n";
   exit;
}

my @args = @ARGV;
if ($args[0] eq "-test") {
   print "Test Mode...\n\n";
   print "Tassel 4 Comment: " . $versionComment4 . "\n";
   print "Tassel 4 Tag: " . $versionTag4 . "\n";
   system "cd ../tassel4.0_standalone; git commit -a -m'$versionComment4' --dry-run";
}
elsif ($args[0] eq "-commit") {
   print "Commit Mode...\n\n";
   print "Tassel 4 Comment: " . $versionComment4 . "\n";
   print "Tassel 4 Tag: " . $versionTag4 . "\n";
   system "cd ../tassel4.0_standalone; git commit -a -m'$versionComment4'";
   system "cd ../tassel4.0_standalone; git tag -m'$versionComment4' $versionTag4";
   system "cd ../tassel4.0_standalone; git push";
   system "cd ../tassel4.0_standalone; git push --tags";
}
