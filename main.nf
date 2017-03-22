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

ref_genome = Channel.fromPath(params.ref_genome).map{ file -> [file.baseName, file] }

process FilterRef {
  tag params.tag
  publishDir "outputs/stages/lastdb"

  cpus 12
  memory { 48.GB }
  time { 6.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  file prefix, fasta from ref_genome

  output:
  file '*.fasta' into filtered_fasta

  when: 

  """
  /bin/zcat $fasta | seqkit seq -m $params.minseqlen > ${prefix}.filtered.fasta 
  """
}

if(params.minseqlen > 0){
  filtered_fasta.into { ref_genome; ref_genome2}
  } else {
  ref_genome.into{ ref_genome; ref_genome2 }
  }

process LastDB {
  tag params.tag
  publishDir "outputs/stages/lastdb"

  cpus 8
  memory { 32.GB }
  time { 6.h }
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
  tag params.tag
  publishDir "outputs/stages/shuffled"

   queue 'dpetrov,normal,hns,owners'
  cpus 1
  memory { 8.GB }
  time { 2.h }
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

shuffled_query.splitFasta( by: 1, file: true, elem: 1).set{split_query}

split_query.into{split_query_view; split_query}

split_query_view.view()

process LastAlign {
  tag params.tag
  publishDir "outputs/stages/alignments"

  cpus 4
  memory { task.attempt == 1 ? 32.GB: 64.GB }
  time { task.attempt == 1 ? 1.d: 2.d }
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
  tag params.tag
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 4.GB }
  time { 1.h}
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

process MakeSamBam {
  tag params.tag
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 4.GB }
  time { 2.h }
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 5
  maxErrors '-1'

  input:
  set qID, file(maf_file) from to_sam
  file ref_fasta from ref_genome2
  output:
  set qID, file("${qID}.sam"), file("${qID}.bam"), file("${qID}.bam.bai") into aligned_sam_bam

  """
  set -e
  /usr/local/bin/maf-convert -n sam ${maf_file} > temp.sam
  /usr/local/bin/samtools faidx ${ref_fasta}
  /usr/local/bin/samtools view -t ${ref_fasta}.fai temp.sam > ${qID}.sam
  /usr/local/bin/samtools view -b -t ${ref_fasta}.fai -T ${ref_fasta} ${qID}.sam > temp.bam
  /usr/local/bin/samtools sort temp.bam > ${qID}.bam
  /usr/local/bin/samtools index ${qID}.bam
  rm temp.sam temp.bam
  """
}

process MakeBlasttab {
  tag params.tag
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 4.GB }
  time { 2.h }
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
  tag params.tag
  publishDir "outputs", mode: 'copy'

  cpus 1
  memory { 4.GB }
  time { 2.h }
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
