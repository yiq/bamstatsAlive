bamstatsAlive
=============

Streaming Bam Stats utility based on Bamtools Stats

It reads a bamfile given on the commandline, or if none is given the stdin, and
output statistics per a given amount of read processed. The default update rate
is 100 reads, but can be changed with the -u parameter

Usage
=====

```
bamstatsalive [options] [bam-file]

Options:
  -u	updateRate [default=100]		The number of reads bamstatsalive needs to process before producing another statistics update
  -f	firstUpdateRate [default=0]		The number of reads bamstatsalive needs to process before producing the first statistics update. Useful to increase app responsiveness
  -r	regionJson	                    A json string describing the sampled regions, needed for coverage histogram. Format: {["chr":"1", "start": 100, "end": 200}, ...]}
  -k	coverageSkipFactor [default=10]	Only 1 in every skipFactor region is used to update coverage histogram, for performance reason.

If no bam-file is specified, input is then read from stdin
```
