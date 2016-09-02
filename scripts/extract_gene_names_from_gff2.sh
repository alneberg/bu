#!/usr/bin/env bash

# Finds the first fasta occurence in the gff2 and only uses lines before that one.
# Finds the values for 'ID=' and 'product=' in the gff product description column

head -n $(( `grep -n '^>' $1 | head -n 1 | cut -f 1 -d ':'` - 2)) $1 | grep -v '^#' | cut -f 9 | sed -E 's/ID=([^;]+).*product=([^;]+).*/\1,\2/g'
