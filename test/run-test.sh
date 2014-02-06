#!/bin/sh

./$1

if [ $? -eq 0 ]; then echo "$1: [32mPASSED[00m"; else echo "$1: [31mFAILED[00m"; fi

exit 0
