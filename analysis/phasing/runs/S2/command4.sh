
docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:dev \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --allele_aware \
--hmmcopyqc SHAH_H003831_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H003831_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts SHAH_H003831_T01_01_DLP01_allele_count.csv.gz \
--ncores 5 \
--rephase_report confidence_4.csv \
--qcplot qcplot_4.png \
--heatmap heatmap_4.png \
--heatmapraw heatmapraw_4.png \
--csvfile hscn_4.csv.gz \
--qccsvfile qc_4.csv \
--Rdatafile signals_4.Rdata \
--mincells 2 \
--sex FEMALE
