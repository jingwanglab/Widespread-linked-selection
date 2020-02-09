#!/usr/bin/awk -f

BEGIN {
	OFS="\t";
}
{
	if ( $0 !~ /\#/ ) {
		print $1, $2-1, $2;
	}
}
