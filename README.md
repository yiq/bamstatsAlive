bamstatsAlive
=============

Streaming Bam Stats utility based on Bamtools Stats

It reads a bamfile given on the commandline, or if none is given the stdin, and
output statistics per a given amount of read processed. The default update rate
is 1000 reads, but can be changed with the -u parameter
