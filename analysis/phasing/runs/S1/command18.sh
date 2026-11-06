
docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:latest \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --phase_llr_min 5 \
  --phase_min_frac 0.5 \
  --relabel_unsupported_phase \
  --rephase_report confidence_18.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_18.png \
--heatmap heatmap_18.png \
--heatmapraw heatmapraw_18.png \
--csvfile hscn_18.csv.gz \
--qccsvfile qc_18.csv \
--Rdatafile signals_18.Rdata \
--mincells 2 \
--sex FEMALE
