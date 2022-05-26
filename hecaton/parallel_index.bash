#!/usr/bin/env bash

module load gnu-parallel/20200722

parallel --block 80M --memfree 20G -a bai.not_finished.txt samtools index
