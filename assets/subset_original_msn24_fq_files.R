# this script is for subsampling the full fastq files to 
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(ShortRead)

set.seed("23423")

# WT mock samples
path_wt_mock_rep1 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 1/JL-AS-24S-TF-11_S11_L002_R1_001.fastq.gz")
path_wt_mock_rep2 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 2/JL-AS-24S-TF-15_S15_L002_R1_001.fastq.gz")
path_wt_mock_rep3 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep3/WT_M_3_JL-AS-72S-pl1-TF-31_S7_L001_R1_001.fastq.gz")
path_wt_mock_rep4 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep4/WT_M_4_JL-AS-72S-pl1-TF-41_S17_L001_R1_001.fastq.gz")

# WT etoh samples
path_wt_etoh_rep1 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 1/JL-AS-24S-TF-08_S8_L002_R1_001.fastq.gz")
path_wt_etoh_rep2 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 2/JL-AS-24S-TF-22_S22_L002_R1_001.fastq.gz")
path_wt_etoh_rep3 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep3/WT_E_3_JL-AS-72S-pl1-TF-25_S1_L001_R1_001.fastq.gz")
path_wt_etoh_rep4 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep4/WT_E_4_JL-AS-72S-pl1-TF-39_S15_L001_R1_001.fastq.gz")

# msn24dd mock samples
path_msn24_mock_rep1 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 1/JL-AS-24S-TF-03_S3_L002_R1_001.fastq.gz")
path_msn24_mock_rep2 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 2/JL-AS-24S-TF-16_S16_L002_R1_001.fastq.gz")
path_msn24_mock_rep3 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep3/msn24_M_3_JL-AS-72S-pl1-TF-29_S5_L001_R1_001.fastq.gz")
path_msn24_mock_rep4 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep4/msn24_M_4_JL-AS-72S-pl1-TF-48_S24_L001_R1_001.fastq.gz")

# msn24dd etoh Samples
path_msn24_etoh_rep1 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 1/JL-AS-24S-TF-02_S2_L002_R1_001.fastq.gz")
path_msn24_etoh_rep2 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep 2/JL-AS-24S-TF-24_S24_L002_R1_001.fastq.gz")
path_msn24_etoh_rep3 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep3/msn24_E_3_JL-AS-72S-pl1-TF-32_S8_L001_R1_001.fastq.gz")
path_msn24_etoh_rep4 <- file.path("~/Library/CloudStorage/Box-Box/Wild Yeast EtOH and NaCl RNA-Seq/TF Mutants/YPS606 TF RNA-Seq v2/YPS606 TF RNA-Seq Rep4/msn24_E_4_JL-AS-72S-pl1-TF-47_S23_L001_R1_001.fastq.gz")


mean_reads <- 200000
sd_reads <- 20000

#wt mock
writeFastq(yield(FastqSampler(path_wt_mock_rep1, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_MOCK_REP1.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_mock_rep2, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_MOCK_REP2.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_mock_rep3, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_MOCK_REP3.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_mock_rep4, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_MOCK_REP4.fastq.gz", mode="w")

#wt etoh
writeFastq(yield(FastqSampler(path_wt_etoh_rep1, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_ETOH_REP1.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_etoh_rep2, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_ETOH_REP2.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_etoh_rep3, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_ETOH_REP3.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_wt_etoh_rep4, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_WT_ETOH_REP4.fastq.gz", mode="w")

#msn24 mock
writeFastq(yield(FastqSampler(path_msn24_mock_rep1, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_MOCK_REP1.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_mock_rep2, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_MOCK_REP2.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_mock_rep3, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_MOCK_REP3.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_mock_rep4, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_MOCK_REP4.fastq.gz", mode="w")

#msn24 etoh
writeFastq(yield(FastqSampler(path_msn24_etoh_rep1, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_ETOH_REP1.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_etoh_rep2, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_ETOH_REP2.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_etoh_rep3, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_ETOH_REP3.fastq.gz", mode="w")
writeFastq(yield(FastqSampler(path_msn24_etoh_rep4, n=rnorm(1,mean=mean_reads,sd=sd_reads))), file="~/Documents/GitHub/GenomicDataAnalysis_Fa23/data/ethanol_stress/fastq/YPS606_MSN24_ETOH_REP4.fastq.gz", mode="w")
