#!/usr/bin/env nextflow

/* Aligns genome 1 to genome 2 using LAST, calculates statistics. See config file for parameters.
  Started October 2016

  @Authors
  Ryan Taylor <ryan@ryantaylor.net>

  To run:
  $ nextflow run align_genomes.nf -c <nextflow.config>
*/


queries = file(params.query)
query_fasta = Channel
  .from(queries.readLines())
  .map {line ->
    lsplit      = line.split()
    queryID     = lsplit[0]
    fastaFile   = file(lsplit[1])
    [ queryID, fastaFile ]
}
query_fasta.into{qf_view; query_fasta}
qf_view.view()

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
  /bin/zcat ${fasta} | /usr/local/bin/lastdb -v -P${task.cpus} ${params.lastdb_options} refdb 
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
  /usr/local/bin/last-train -P${task.cpus} ldb ${fasta1} > genome1.par
  """
}
*/

process ShuffleFasta {
  publishDir "outputs/stages/shuffled"

   queue 'dpetrov,normal,hns,owners'
  cpus 2
  memory { 4.GB * task.cpus }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set qID, file(fasta) from query_fasta

  output:
  set qID, file("*.shuf.fasta") into shuffled_query

  """
  /usr/local/bin/seqkit shuffle ${fasta} > ${qID}.shuf.fasta
  """
}

shuffled_query.splitFasta( by: 100, file: true, elem: 1).set{split_query}

split_query.into{split_query_view; split_query}

split_query_view.view()

process LastAlign {
  publishDir "outputs/stages/alignments"

  cpus 4
  memory { 32.GB * task.attempt }
  time { 2.d }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 4
  maxErrors '-1'
 
  input:
  set qID, file(fasta_i) from split_query
  file db_files from database_files.first()
 
  output:
  set qID, file("*.maf") into aligned_mafs
 
  """
  /usr/local/bin/lastal -v -P${task.cpus} ${params.lastal_options} refdb ${fasta_i} | /usr/local/bin/last-split -v -m1 > ${fasta_i}.maf
  """
}

aligned_mafs.groupTuple().into{aligned_mafs_view; aligned_mafs}

aligned_mafs_view.view()

process MergeMAF {
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 32.GB * task.cpus }
  time { 2.d * task.attempt }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set qID, file(mafs) from aligned_mafs

  output:
  set qID, file("*.maf") into aligned_maf

  script:
  input_mafs = mafs.collect{"$it"}.join(' ')

  """
  /bin/cat ${input_mafs} > ${qID}.maf
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
  set qID, file(maf_file) from to_sam

  output:
  set qID, file("*.sam") into aligned_sam

  """
  /usr/local/bin/maf-convert -n sam ${maf_file} > ${qID}.sam
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
  set qID, file(sam_file) from aligned_sam
  file ref_file from ref_genome2

  output:
  set qID, file("*.bam") into aligned_bam

  """
  /usr/local/bin/samtools view -bT ${ref_file} ${sam_file} > ${qID}.bam
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
  set qID, file(bam_file) from aligned_bam

  output:
  set qID, file("*.sorted.bam") into sorted_sam

  """
  /usr/local/bin/samtools sort -T samtools_tmp ${bam_file} > ${qID}.sorted.bam
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
  set qID, file(maf_file) from to_blasttab

  output:
  set qID, file("*.blasttab") into aligned_blasttab

  """
  /usr/local/bin/maf-convert -n blasttab ${maf_file} > ${qID}.blasttab
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
  set qID, file(maf_file) from to_tab

  output:
  set qID, file("*.tab") into aligned_tab

  """
  /usr/local/bin/maf-convert -n tab ${maf_file} > ${qID}.tab
  """
}
