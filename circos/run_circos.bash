#!/usr/bin/env bash

# singularity run -B /zfs ~/singularity_containers/circos_latest.sif -conf etc/circos_sap.conf
singularity run -B /zfs ~/singularity_containers/circos_latest.sif -conf etc/circos_variants.conf
