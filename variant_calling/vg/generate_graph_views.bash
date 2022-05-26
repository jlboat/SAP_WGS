#!/usr/bin/env bash

while read p
do
    vg view ${p}.vg > ${p}.gfa
done < <(ls -1 chunk_*vg | cut -f 1 -d '.')


# chunk_0_Chr09_57038604_57141201.vg > chunk_0_Chr09_57038604_57141201.gfa
