#!/usr/bin/env nextflow

/* Aligns genome 1 to genome 2 using LAST, calculates statistics. See config file for parameters.
  Started October 2016

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

  To run:
  $ nextflow run align_genomes.nf -c <nextflow.config>
*/

/*process TestPipeline {
  publishDir "outputs/tests"

   queue 'dpetrov,normal,hns,owners'
  cpus 2
  memory { 4.GB * task.cpus }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file fasta1 from genome1
  file fasta2 from genome2

  output:
  file "sample1.fasta" into testgenome1
  file "sample2.fasta" into testgenome2

  """
  parallel -j${task.cpus} --link 'zcat {1} | seqkit sample -n 1000 -o {2}' ::: ${fasta1} ${fasta2} ::: sample1.fasta sample2.fasta
  """
}
*/

Channel.fromPath(params.query_fastas).set{query_fastas}
Channel.fromPath(params.ref_genome).set{ref_genome}
ref_genome.into{ref_genome; ref_genome2}

process LastDB {
  publishDir "outputs/stages/lastdb"

  cpus 16
  memory { 4.GB * task.cpus }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file fasta from ref_genome

  output:
  file 'refdb*' into database_files

  """
  zcat ${fasta} | lastdb -v -P${task.cpus} ${params.lastdb_options} refdb 

  """
}

/*process LastTrain {
  publishDir "outputs/stages/training"

  cpus 32
  memory { 4.GB * task.cpus }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file fasta1 from testgenome1
  file(database_files)

  output:
  file 'genome1.par' into database_params

  """
  last-train -P${task.cpus} ldb ${fasta1} > genome1.par
  """
}
*/

Channel.fromPath(params.query_fastas).splitFasta( by: 4, file: true).set{split_query}

process LastAlign {
  publishDir "outputs/stages/alignments"

  cpus 4
  memory { 32.GB * task.attempt }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'
 
  input:
  file fasta_i from split_query
  file db_files from database_files.first()
 
  output:
  file "${fasta_i}.maf" into align_mafs
 
  """
  lastal -v -P${task.cpus} ${params.lastal_options} refdb ${fasta_i} | last-split -v -m1 > ${fasta_i}.maf
  """
}

align_mafs = align_mafs.toList()

process MergeMAF {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file mafs from align_mafs

  output:
  file "${params.output_prefix}.maf" into aligned_maf

  script:
  input_mafs = mafs.collect{"$it"}.join(' ')

  """
  cat ${input_mafs} > ${params.output_prefix}.maf
  """
}

aligned_maf.into{to_sam; to_blasttab; to_tab}

process MakeSam {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file maf_file from to_sam

  output:
  file "${params.output_prefix}.sam" into aligned_sam

  """
  maf-convert -n sam ${maf_file} > ${params.output_prefix}.sam
  """
}

process ConvertSamToBam {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file sam_file from aligned_sam
  file ref_file from ref_genome2

  output:
  file "${params.output_prefix}.sam" into aligned_bam

  """
  samtools view -bT ${ref_file} ${sam_file} > ${params.output_prefix}.sam
  """
}

process SortBam {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file bam_file from aligned_bam

  output:
  file "${params.output_prefix}.sam" into sorted_sam

  """
  samtools sort -T samtools_tmp ${bam_file} > ${params.output_prefix}.sorted.bam
  """
}

process MakeBlasttab {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file maf_file from to_blasttab

  output:
  file "${params.output_prefix}.blasttab" into aligned_blasttab

  """
  maf-convert -n blasttab ${maf_file} > ${params.output_prefix}.blasttab
  """
}

process MakeTab {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file maf_file from to_tab

  output:
  file "${params.output_prefix}.tab" into aligned_tab

  """
  maf-convert -n tab ${maf_file} > ${params.output_prefix}.tab
  """
}








