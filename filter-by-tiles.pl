#!/usr/bin/perl -w

$file = $ARGV[0];

#print $file,"\n";

$count = "80";
if ($#ARGV > 0) {
    $count = $ARGV[1];
}
open FILE,$file or die "Can't find $file\n";

while (<FILE>) {
    @fields = split /\s+/;
    $size = @fields;
    if ($size > 1 && $size < 12 && $fields[5] > $count) {
	print $_;
    }	
}
