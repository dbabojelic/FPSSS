#!/usr/bin/perl
#
#

CON: while (<>) {
    chomp;
    if (/^>/) {
        print;
        print "\n";
        next CON;
    }

    $len = length;
    $start = 0;
    $end = 60;
    if ($len < $end) {
        $end = $len;
    }

    while ($end < $len) {
        print substr $_, $start, ($end - $start);
        print "\n";
        $start += 60;
        $end += 60;
        if ($len < $end) {
            $end = $len;
        }
    }
    print substr $_, $start, ($end - $start);
    print "\n";
}
