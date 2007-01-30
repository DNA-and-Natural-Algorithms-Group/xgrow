#!/usr/bin/perl -w

$file = $ARGV[0];

print $file,"\n";

$count = "80";
if ($#ARGV > 0) {
    $count = $ARGV[1];
}
open FILE,$file or die "Can't find $file\n";

while (<FILE>) {
    @fields = split /\s+/;
    if ($fields[5] eq $count) {
	print $_;
    }
}
