
docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --rephase_by_imbalance \
  --relabel_unsupported_phase \
  --rephase_report rephasing_13.csv \
--hmmcopyqc SHAH_H003831_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H003831_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts SHAH_H003831_T01_01_DLP01_allele_count.csv.gz \
--ncores 5 \
--qcplot qcplot_1.png \
--heatmap heatmap_1.png \
--heatmapraw heatmapraw_1.png \
--csvfile hscn_1.csv.gz \
--qccsvfile qc_1.csv \
--Rdatafile signals_1.Rdata \
--mincells 2 \
--sex FEMALE

