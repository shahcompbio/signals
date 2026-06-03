sim_data_bb <- simulate_data_cohort(
  clone_num = c(20, 25, 25, 10),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("17" = c(3, 1), "8" = c(6, 2)),
    list("1" = c(2, 2), "9" = c(4, 1))
  ), # opposite LOH on chr 1
  loherror = 0.01,
  coverage = 100,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

dfannot <- data.frame(cell_id = unique(sim_data_bb$ascn$cell_id)) 
dfannot$Event <- NULL
dfannot$Event[1:20] <- "1"
dfannot$Event[21:45] <- "2"
dfannot$Event[46:70] <- "3"
dfannot$Event[71:80] <- "4"
dfannot$Clone <- "0"
dfannot$Other <- NULL
dfannot$Other[1:40] <- "Other1"
dfannot$Other[41:80] <- "Other2"

hm1 <- plotHeatmap(sim_data_bb$ascn, 
            annotations = dfannot, 
            tree = NULL, 
            reorderclusters = TRUE, 
            plottree = FALSE)

cl <- umap_clustering(sim_data_bb$ascn, 
                     field = "state",
                     minPts = 5,
                     umapmetric = "euclidean")
dfannot2 <- dplyr::left_join(dfannot, 
                            cl$clustering[c("cell_id", "clone_id")], 
                            by = "cell_id")
hm2 <- plotHeatmap(sim_data_bb$ascn, 
            annotations = dfannot2, 
            clusters = cl$clustering,
            tree = cl$tree, 
            reorderclusters = TRUE, 
            plottree = TRUE)

dfannot_continuous <- data.frame(cell_id = unique(sim_data_bb$ascn$cell_id))
dfannot_continuous$RedScore <- seq(0, 1, length.out = nrow(dfannot_continuous))
dfannot_continuous$BlueScore <- seq(10, 89, length.out = nrow(dfannot_continuous))
dfannot_continuous$PurpleScore <- seq(100, 179, length.out = nrow(dfannot_continuous))
dfannot_continuous$OrangeScore <- seq(1000, 1079, length.out = nrow(dfannot_continuous))

dfannot_threshold <- data.frame(cell_id = unique(sim_data_bb$ascn$cell_id))
dfannot_threshold$Score <- rep(seq_len(10), length.out = nrow(dfannot_threshold))

dfannot_manual <- data.frame(
  cell_id = unique(sim_data_bb$ascn$cell_id),
  Group = dfannot$Event,
  Score = dfannot_continuous$RedScore
)

manual_annotation_colours <- list(
  Group = c(
    `1` = "#1B9E77",
    `2` = "#D95F02",
    `3` = "#7570B3",
    `4` = "#E7298A"
  ),
  Score = circlize::colorRamp2(c(0, 1), c("#F5F4F0", "#7A3E9D"))
)

hm_continuous <- plotHeatmap(
  sim_data_bb$ascn,
  annotations = dfannot_continuous,
  tree = NULL,
  reorderclusters = TRUE,
  plottree = FALSE
)

hm_threshold_default <- plotHeatmap(
  sim_data_bb$ascn,
  annotations = dfannot_threshold,
  tree = NULL,
  reorderclusters = TRUE,
  plottree = FALSE
)

hm_threshold <- plotHeatmap(
  sim_data_bb$ascn,
  annotations = dfannot_threshold,
  annotation_continuous_threshold = 9,
  tree = NULL,
  reorderclusters = TRUE,
  plottree = FALSE
)

hm_manual_colours <- plotHeatmap(
  sim_data_bb$ascn,
  annotations = dfannot_manual,
  annotation_colours = manual_annotation_colours,
  tree = NULL,
  reorderclusters = TRUE,
  plottree = FALSE
)

hm_meaniqr <- plotHeatmap(
  sim_data_bb$ascn,
  clusters = cl$clustering,
  plotcol = "copy",
  plotmeaniqr = TRUE,
  reorderclusters = TRUE,
  plottree = FALSE
)

hm_meaniqr_from_copy <- plotHeatmap(
  sim_data_bb$ascn,
  clusters = cl$clustering,
  plotcol = "state",
  plotmeaniqr = TRUE,
  meaniqr_plotcol = "copy",
  reorderclusters = TRUE,
  plottree = FALSE
)

copynumber_copy <- signals:::createCNmatrix(sim_data_bb$ascn, field = "copy")
copynumber_copy_formatted <- signals:::format_copynumber(
  copynumber_copy,
  ordered_cell_ids = unique(sim_data_bb$ascn$cell_id),
  plotcol = "copy"
)

test_that("Test returns plot object", {
  expect_true(typeof(hm1) == "S4")
  expect_true(typeof(hm2) == "S4")
  expect_true(typeof(hm_continuous) == "S4")
  expect_true(typeof(hm_threshold_default) == "S4")
  expect_true(typeof(hm_threshold) == "S4")
  expect_true(typeof(hm_manual_colours) == "S4")
  expect_true(typeof(hm_meaniqr) == "S4")
  expect_true(typeof(hm_meaniqr_from_copy) == "S4")
})

test_that("Continuous annotations use pale to accent gradients in order", {
  ha <- signals:::make_left_annot_generic(dfannot_continuous)

  expect_identical(
    ha@anno_list$RedScore@color_mapping@col_fun(c(0, 1)),
    c("#F5F4F0FF", "#C95D63FF")
  )
  expect_identical(
    ha@anno_list$BlueScore@color_mapping@col_fun(c(10, 89)),
    c("#F5F4F0FF", "#5B84B1FF")
  )
  expect_identical(
    ha@anno_list$PurpleScore@color_mapping@col_fun(c(100, 179)),
    c("#F5F4F0FF", "#8C6BB1FF")
  )
  expect_identical(
    ha@anno_list$OrangeScore@color_mapping@col_fun(c(1000, 1079)),
    c("#F5F4F0FF", "#D99A4EFF")
  )
})

test_that("Continuous annotation threshold is configurable", {
  ha_default <- signals:::make_left_annot_generic(dfannot_threshold)
  ha_custom <- signals:::make_left_annot_generic(
    dfannot_threshold,
    continuous_threshold = 9
  )

  expect_identical(ha_default@anno_list$Score@color_mapping@type, "discrete")
  expect_identical(ha_custom@anno_list$Score@color_mapping@type, "continuous")
})

test_that("Annotation colour overrides support discrete and continuous columns", {
  ha_manual <- signals:::make_left_annot_generic(
    dfannot_manual,
    annotation_colours = manual_annotation_colours
  )

  expect_identical(
    ha_manual@anno_list$Group@color_mapping@full_col,
    c(
      `1` = "#1B9E77FF",
      `2` = "#D95F02FF",
      `3` = "#7570B3FF",
      `4` = "#E7298AFF"
    )
  )
  expect_identical(
    ha_manual@anno_list$Score@color_mapping@col_fun(c(0, 1)),
    c("#F5F4F0FF", "#7A3E9DFF")
  )
})

test_that("Annotation colour overrides validate discrete mappings", {
  expect_error(
    signals:::make_left_annot_generic(
      dfannot_manual,
      annotation_colours = list(Group = c(`1` = "#1B9E77"))
    ),
    "does not cover all values"
  )
})

test_that("Mean plus IQR summary annotation renders", {
  summary_annot <- signals:::make_summary_annotations(
    copynumber_copy_formatted,
    plotcol = "copy",
    plotmeaniqr = TRUE
  )

  expect_true(typeof(summary_annot) == "S4")
  expect_true("mean_iqr_cn" %in% names(summary_annot@anno_list))
})

test_that("Mean plus IQR track can use a different source column", {
  summary_annot <- signals:::make_summary_annotations(
    copynumber = hm1@matrix,
    meaniqr_copynumber = copynumber_copy_formatted,
    plotcol = "state",
    meaniqr_plotcol = "copy",
    plotmeaniqr = TRUE
  )

  expect_true(typeof(summary_annot) == "S4")
  expect_true("mean_iqr_cn" %in% names(summary_annot@anno_list))
  expect_true("mean_iqr_cn" %in% names(hm_meaniqr_from_copy@top_annotation@anno_list))
})

test_that("Missing annotation cells are removed before reordering", {
  dfannot_missing <- dfannot[-1, , drop = FALSE]

  expect_warning(
    hm_missing <- plotHeatmap(
      sim_data_bb$ascn,
      annotations = dfannot_missing,
      tree = cl$tree,
      plottree = FALSE
    ),
    "removing non-overlapping cells"
  )

  expect_true(typeof(hm_missing) == "S4")
})

# Gene annotation tests
test_that("Gene annotations render correctly", {
  # Test with valid genes on chromosomes in the data
  hm_genes <- plotHeatmap(sim_data_bb$ascn,
                          gene_annotations = c("TP53", "MYC"),
                          plottree = FALSE,
                          reorderclusters = TRUE)
  expect_true(typeof(hm_genes) == "S4")
})

test_that("Gene annotations work with frequency plot", {
  hm_genes_freq <- plotHeatmap(sim_data_bb$ascn,
                               gene_annotations = c("TP53"),
                               plotfrequency = TRUE,
                               plottree = FALSE,
                               reorderclusters = TRUE)
  expect_true(typeof(hm_genes_freq) == "S4")
})

test_that("Invalid gene names produce warning", {
  expect_warning(
    plotHeatmap(sim_data_bb$ascn,
                gene_annotations = c("FAKE_GENE_NAME"),
                plottree = FALSE,
                reorderclusters = TRUE),
    "Gene\\(s\\) not found"
  )
})

test_that("Gene annotation parameters are respected", {
  hm_custom <- plotHeatmap(sim_data_bb$ascn,
                           gene_annotations = c("TP53", "MYC"),
                           gene_annotation_fontsize = 8,
                           gene_link_height = 10,
                           gene_label_sep = " | ",
                           plottree = FALSE,
                           reorderclusters = TRUE)
  expect_true(typeof(hm_custom) == "S4")
})
