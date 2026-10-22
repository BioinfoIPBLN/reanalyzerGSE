#!/bin/bash
exec awk -F '\t' '
function flush(   k) {
	if (n == 0) return
	if (single || (r1 && r2)) { for (k = 1; k <= n; k++) print buf[k]; kept++ }
	else dropped++
	n = 0; r1 = 0; r2 = 0; single = 0
}
/^@/ { print; next }
$1 != q { flush(); q = $1 }
{
	buf[++n] = $0
	f = $2 + 0
	if (f % 2 == 0) single = 1
	else if (int(f / 256) % 2 == 0 && int(f / 2048) % 2 == 0) {
		if (int(f / 64) % 2 == 1) r1 = 1
		if (int(f / 128) % 2 == 1) r2 = 1
	}
}
END {
	flush()
	printf "sam_keep_pairs: %d reads or pairs kept, %d pairs dropped because only one mate passed the BAM filters\n", kept, dropped > "/dev/stderr"
}'
