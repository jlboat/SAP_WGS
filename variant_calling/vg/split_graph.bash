#!/usr/bin/env bash

# vg index -x wg_pruned.xg $(for i in $(seq 1 9;); do echo Chr0${i}.pruned.vg; done) Chr10.pruned.vg
vg chunk -x wg.xg --input-bed sap_genes.bed -c 1
