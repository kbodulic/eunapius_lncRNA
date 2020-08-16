use strict;
use warnings;

my @arr;

while (<>) {
	    chomp;
	        push @arr, $_ if length;
		    last if eof;
	    }

	    while (<>) {
		        print /^>/ ? shift(@arr) . "\n" : $_;
		}
