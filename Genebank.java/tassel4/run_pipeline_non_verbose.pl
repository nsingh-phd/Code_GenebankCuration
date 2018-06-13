#!/usr/bin/perl -w

use strict;
my $top = '.';
my $libdir = "$top/lib";
opendir (DIR, "$libdir") || die "Could not open $libdir\n";
my @list = readdir(DIR);

my @fl = ();
foreach my $fn(@list){
   if ("$fn" =~ m/\.jar$/){
      push(@fl, "$libdir\/$fn");
   }
}
push(@fl, "$top/dist/sTASSEL.jar");
my $CP = join(":", @fl);
my @args = @ARGV;

system "java -classpath '$CP' -Xms512m -Xmx4096m net.maizegenetics.pipeline.TasselPipeline @args";
