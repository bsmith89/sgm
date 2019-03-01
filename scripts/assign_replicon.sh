#!/usr/bin/env bash


join -j 2 \
    <(cat $3 \
          | awk -v min_ident=$1 -v min_len=$2 '$3 > min_ident && $4 > min_len {print $1,$2}' \
          | sort -k2,2) \
    <(cat $4 \
          | sed '1,1d' \
          | cut -f1,3 \
          | sort -k2,2) \
    | sort | uniq | awk -v OFS='\t' '{print $2,$3}'
