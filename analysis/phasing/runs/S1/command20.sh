
docker run --platform linux/x86_64 -w `pwd` -v /Users/mcphera1/Projects/signals:/Users/mcphera1/Projects/signals signals:dev \
Rscript run_signals.R \
  --seed 42 \
  --viterbiver cpp \
  --rephase_unit arm \
  --joint_phasing \
  --joint_iters 4 \
  --phase_llr_min 5 \
  --phase_min_frac 0.5 \
  --allele_aware \
  --rephase_report confidence_20.csv \
--hmmcopyqc SHAH_H002216_T01_01_DLP01_hmmcopy_metrics.csv.gz \
--hmmcopyreads SHAH_H002216_T01_01_DLP01_hmmcopy_reads.csv.gz \
--allelecounts haplotype_calling_count.csv.gz \
--ncores 5 \
--qcplot qcplot_20.png \
--heatmap heatmap_20.png \
--heatmapraw heatmapraw_20.png \
--csvfile hscn_20.csv.gz \
--qccsvfile qc_20.csv \
--Rdatafile signals_20.Rdata \
--mincells 2 \
--sex FEMALE
