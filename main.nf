#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DAREVSKIA } from './workflows/darevskia.nf'

workflow {
    DAREVSKIA()
}
