#!/usr/bin/perl
#

$cnt = 0;
while (<>) {
    chomp;
    ($id) = split /\s/;
    if ($id eq $lastId) {
        $cnt++;
    }
    else {
        $cnt = 0;
    }
    if ($cnt < 5) {
        print;
        print "\n";
    }
    $lastId = $id;
}
