#!/bin/sh

export X3DNA=/home/cosi/Bio/X3DNA
export PATH=$X3DNA/bin:$PATH

find_pair $1 $2 > /dev/null 2>&1

