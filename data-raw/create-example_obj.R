library(TraceQC)
library(fastqcr)
input_file <- system.file("extdata", "test_data",
                          "fastq", "example.fastq.gz", package="TraceQC")
ref_file <- system.file("extdata", "test_data", "ref",
                        "ref.txt", package="TraceQC")
qc_dir <- tempdir()
fastqc(system.file("extdata", "test_data",
                   "fastq", package = "TraceQC"),
       qc.dir=qc_dir)

input_qc_path <- get_qcpath(input_file, qc_dir)

obj <- TraceQC(input_file = input_file,
               ref_file = ref_file,
               fastqc_file = input_qc_path,
               ncores=1)


example_obj <- obj
save(example_obj, file = "data/example_obj.rda")
