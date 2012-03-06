#!/bin/bash -x
for xdtype in sr dr sc dc ; do
	 for ydtype in sr dr sc dc ; do
		   ./tprod_testcases $xdtype 11x29x51 rand '-1 1 3' \
				$ydtype 11x29x51 rand '-1 2 3' #> /dev/null
		   ./tprod_testcases $xdtype 11x29x51 rand '-1 2 3' \
				$ydtype 11x29x51 rand '-1 1 3' #> /dev/null
		   ./tprod_testcases $xdtype 11x29x51 rand '1 -2 3' \
				$ydtype 11x29x51 rand '2 -2 3' #> /dev/null
		   ./tprod_testcases $xdtype 11x29x51 rand '2 -2 3' \
				$ydtype 11x29x51 rand '1 -2 3' #> /dev/null
	 done
done
