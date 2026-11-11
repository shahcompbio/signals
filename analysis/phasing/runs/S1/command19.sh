
docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --phase_llr_min 5 \
  --phase_min_frac 0.5 \
  --rephase_report confidence_19.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_19.png \
--heatmap heatmap_19.png \
--heatmapraw heatmapraw_19.png \
--csvfile hscn_19.csv.gz \
--qccsvfile qc_19.csv \
--Rdatafile signals_19.Rdata \
--mincells 2 \
--sex FEMALE
