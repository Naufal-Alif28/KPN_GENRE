#!/usr/bin/env bash

mkdir -p "cache/sra"

while IFS= read -r sra; do
    prefetch ${sra} --output-file "cache/sra/${sra}.sra"
done < "sra_list.txt"