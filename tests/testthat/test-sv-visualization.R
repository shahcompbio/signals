library(testthat)
library(ggplot2)
library(dplyr)

# Create mock CNbins data
set.seed(123)
mock_CNbins <- data.frame(
  cell_id = rep("test_cell", 100),
  chr = rep(c("11", "6"), each = 50),
  start = c(seq(1, 50e6, by = 1e6), seq(1, 50e6, by = 1e6)),
  end = c(seq(1e6, 50e6, by = 1e6), seq(1e6, 50e6, by = 1e6)),
  state = sample(0:4, 100, replace = TRUE),
  copy = rnorm(100, mean = 2, sd = 0.5)
)

# Create mock SV data with all orientation types
mock_SV <- data.frame(
  chromosome_1 = c("11", "11", "11", "11", "11"),
  chromosome_2 = c("11", "11", "11", "11", "6"),
  position_1 = c(10e6, 20e6, 30e6, 40e6, 50e6),
  position_2 = c(15e6, 25e6, 35e6, 45e6, 10e6),
  strand_1 = c("+", "-", "+", "-", "+"),
  strand_2 = c("-", "+", "+", "-", "-"),
  read_count = c(100, 150, 200, 120, 180),
  rearrangement_type = c("deletion", "duplication", "inversion", "deletion", "translocation"),
  type = c("DEL", "DUP", "INV", "DEL", "TRA")
)

test_that("classify_sv_orientation returns correct categories", {
  # Test intra-chromosomal orientations
  sv_intra <- data.frame(
    chromosome_1 = c("11", "11", "11", "11"),
    chromosome_2 = c("11", "11", "11", "11"),
    strand_1 = c("+", "-", "+", "-"),
    strand_2 = c("-", "+", "+", "-")
  )
  
  orientations <- classify_sv_orientation(sv_intra)
  expect_equal(orientations, c("+-", "-+", "++", "--"))
  
  # Test inter-chromosomal (translocation)
  sv_inter <- data.frame(
    chromosome_1 = c("11", "6"),
    chromosome_2 = c("6", "11"),
    strand_1 = c("+", "-"),
    strand_2 = c("-", "+")
  )
  
  orientations <- classify_sv_orientation(sv_inter)
  expect_equal(orientations, c("Translocation", "Translocation"))
  
  # Test empty input
  expect_equal(classify_sv_orientation(data.frame()), character(0))
  expect_equal(classify_sv_orientation(NULL), character(0))
})

test_that("prepare_sv_points creates correct structure", {
  binsize <- 1e6
  bins <- data.frame(
    chr = rep(c("11", "6"), each = 50),
    start = c(seq(1, 50e6, by = 1e6), seq(1, 50e6, by = 1e6)),
    idx = 1:100
  )
  
  sv_points <- prepare_sv_points(mock_SV, bins, binsize)
  
  expect_true(nrow(sv_points) > 0)
  expect_true(all(c("idx", "read_count", "orientation", "breakpoint_id") %in% names(sv_points)))
  expect_true(is.numeric(sv_points$idx))
  expect_true(is.numeric(sv_points$read_count))
  expect_true(is.character(sv_points$orientation))
  
  # Test with empty SV
  empty_sv <- data.frame(
    chromosome_1 = character(0),
    chromosome_2 = character(0),
    position_1 = numeric(0),
    position_2 = numeric(0),
    strand_1 = character(0),
    strand_2 = character(0),
    read_count = numeric(0)
  )
  empty_points <- prepare_sv_points(empty_sv, bins, binsize)
  expect_equal(nrow(empty_points), 0)
})

test_that("generate_sv_arcs creates valid paths", {
  sv_with_idx <- data.frame(
    idx_1 = c(10, 20, 30),
    idx_2 = c(15, 25, 35),
    orientation = c("+-", "-+", "++")
  )
  
  y_start <- c(1, 1.5, 2)
  y_end <- c(1.2, 1.7, 2.2)
  
  arcs <- generate_sv_arcs(
    sv_with_idx,
    y_start = y_start,
    y_end = y_end,
    arc_height_factor = 0.5,
    n_points = 10
  )
  
  expect_true(nrow(arcs) > 0)
  expect_true(all(c("idx", "y", "arc_id", "orientation") %in% names(arcs)))
  expect_true(is.numeric(arcs$idx))
  expect_true(is.numeric(arcs$y))
  expect_true(all(arcs$y >= 0))  # Arcs should be above x-axis
  
  # Test with same start/end (should be filtered out)
  sv_same <- data.frame(
    idx_1 = c(10, 20),
    idx_2 = c(10, 20),
    orientation = c("+-", "-+")
  )
  arcs_same <- generate_sv_arcs(
    sv_same,
    y_start = c(1, 1),
    y_end = c(1.2, 1.2),
    arc_height_factor = 0.5
  )
  expect_equal(nrow(arcs_same), 0)
  
  # Test with empty input
  empty_arcs <- generate_sv_arcs(
    data.frame(),
    y_start = numeric(0),
    y_end = numeric(0),
    arc_height_factor = 0.5
  )
  expect_equal(nrow(empty_arcs), 0)
})

test_that("plotCNprofile with sv_style = lines_and_arcs works", {
  # Test basic functionality
  p <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "lines_and_arcs",
    chrfilt = c("11", "6")
  )
  
  expect_true(inherits(p, "ggplot"))
  
  # Test with both styles
  p_both <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "both",
    chrfilt = c("11", "6")
  )
  
  expect_true(inherits(p_both, "ggplot"))
  
  # Test backward compatibility (curves style)
  p_curves <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "curves",
    chrfilt = c("11", "6")
  )
  
  expect_true(inherits(p_curves, "ggplot"))
})

test_that("plotCNprofile squashy axis with SV read axis works", {
  p_squashy <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "lines_and_arcs",
    chrfilt = c("11", "6"),
    y_axis_trans = "squashy",
    show_sv_read_axis = TRUE
  )

  expect_true(inherits(p_squashy, "ggplot"))
})

test_that("plotCNprofile validation catches invalid sv_style", {
  expect_error(
    plotCNprofile(mock_CNbins, SV = mock_SV, sv_style = "invalid"),
    "sv_style must be one of"
  )
})

test_that("plotCNprofile validation catches missing strand columns", {
  sv_no_strand <- mock_SV %>% dplyr::select(-strand_1, -strand_2)
  
  expect_error(
    plotCNprofile(mock_CNbins, SV = sv_no_strand, sv_style = "lines_and_arcs"),
    "SV data missing required columns"
  )
  
  sv_no_readcount <- mock_SV %>% dplyr::select(-read_count)
  
  expect_error(
    plotCNprofile(mock_CNbins, SV = sv_no_readcount, sv_style = "lines_and_arcs"),
    "SV data missing required columns"
  )
})

test_that("SV_orientation_cols function works", {
  colors <- SV_orientation_cols()
  expect_true(is.character(colors))
  expect_true(length(colors) >= 5)
  expect_true(all(c("+-", "-+", "++", "--", "Translocation") %in% names(colors)))

  # Test subsetting
  subset_colors <- SV_orientation_cols("+-", "-+")
  expect_equal(length(subset_colors), 2)
})

test_that("plotCNprofile with lines_and_arcs and squashy transform works", {
  # Test that squashy + secondary axis produces no errors/warnings
  expect_no_error({
    p <- plotCNprofile(
      mock_CNbins,
      cellid = "test_cell",
      SV = mock_SV,
      sv_style = "lines_and_arcs",
      y_axis_trans = "squashy",
      show_sv_read_axis = TRUE,
      chrfilt = c("11", "6")
    )
  })

  p <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "lines_and_arcs",
    y_axis_trans = "squashy",
    show_sv_read_axis = TRUE,
    chrfilt = c("11", "6")
  )

  expect_true(inherits(p, "ggplot"))

  # Test with identity transform as well (should still work)
  p_identity <- plotCNprofile(
    mock_CNbins,
    cellid = "test_cell",
    SV = mock_SV,
    sv_style = "lines_and_arcs",
    y_axis_trans = "identity",
    show_sv_read_axis = TRUE,
    chrfilt = c("11", "6")
  )

  expect_true(inherits(p_identity, "ggplot"))
})

test_that("get_sv_lines_and_arcs_legend returns a legend grob", {
  leg <- get_sv_lines_and_arcs_legend()

  # Check it returns a grob
  expect_s3_class(leg, "gtable")

  # Check it's not empty
  expect_true(length(leg$grobs) > 0)
})

test_that("get_sv_lines_and_arcs_legend accepts custom parameters", {
  leg <- get_sv_lines_and_arcs_legend(
    legend_title = "Custom Title",
    text_size = 12,
    title_size = 14
  )

  expect_s3_class(leg, "gtable")
})

test_that("get_sv_lines_and_arcs_legend can be combined with plotCNprofile", {
  # Create plot without legend using mock data
  p <- plotCNprofile(mock_CNbins,
                     cellid = "test_cell",
                     SV = mock_SV,
                     sv_style = "lines_and_arcs",
                     legend.position = "none",
                     chrfilt = c("11", "6"))

  # Get legend
  leg <- get_sv_lines_and_arcs_legend()

  # Combine with cowplot
  combined <- cowplot::plot_grid(p, leg, rel_widths = c(1, 0.2))

  expect_s3_class(combined, "ggplot")
})

test_that("get_sv_lines_and_arcs_legend supports horizontal direction", {
  # Test vertical (default)
  leg_vert <- get_sv_lines_and_arcs_legend(direction = "vertical")
  expect_s3_class(leg_vert, "gtable")

  # Test horizontal
  leg_horiz <- get_sv_lines_and_arcs_legend(direction = "horizontal")
  expect_s3_class(leg_horiz, "gtable")

  # Both should create valid legend objects with non-zero dimensions
  # The direction parameter affects internal layout but may not change gtable dimensions
  expect_true(all(dim(leg_vert) > 0))
  expect_true(all(dim(leg_horiz) > 0))
})

test_that("get_sv_lines_and_arcs_legend validates direction parameter", {
  expect_error(
    get_sv_lines_and_arcs_legend(direction = "invalid"),
    "direction must be either 'vertical' or 'horizontal'"
  )
})

test_that("flip_sv_positions handles intra-chromosomal SVs with reversed positions", {
  # Create mock SV data with reversed positions
  mock_SV_reversed <- data.frame(
    chromosome_1 = c("11", "11", "11", "6"),
    chromosome_2 = c("11", "11", "6", "11"),  # Mix of intra and inter-chromosomal
    position_1 = c(20e6, 30e6, 10e6, 5e6),    # Two need flipping
    position_2 = c(10e6, 25e6, 20e6, 15e6),   # position_1 > position_2 for rows 1,2
    strand_1 = c("+", "-", "+", "-"),
    strand_2 = c("-", "+", "-", "+"),
    read_count = c(100, 150, 200, 120),
    rearrangement_type = c("deletion", "duplication", "translocation", "deletion"),
    type = c("DEL", "DUP", "TRA", "DEL")
  )

  # Should give warning
  expect_warning(
    result <- signals:::flip_sv_positions(mock_SV_reversed),
    "Flipped positions and strands for 2 intra-chromosomal SV"
  )

  # Check that positions are now correct for intra-chromosomal SVs
  intra <- result[result$chromosome_1 == result$chromosome_2, ]
  expect_true(all(intra$position_1 < intra$position_2))

  # Check that strands were also flipped for row 1
  expect_equal(result$strand_1[1], "-")  # Was "+"
  expect_equal(result$strand_2[1], "+")  # Was "-"

  # Check that strands were also flipped for row 2
  expect_equal(result$strand_1[2], "+")  # Was "-"
  expect_equal(result$strand_2[2], "-")  # Was "+"

  # Check that inter-chromosomal SV (row 3) was NOT touched
  expect_equal(result$position_1[3], 10e6)
  expect_equal(result$position_2[3], 20e6)
  expect_equal(result$strand_1[3], "+")
  expect_equal(result$strand_2[3], "-")
})

test_that("flip_sv_positions handles already-correct positions", {
  # Create mock SV data with correct positions (position_1 < position_2)
  mock_SV_correct <- data.frame(
    chromosome_1 = c("11", "11"),
    chromosome_2 = c("11", "11"),
    position_1 = c(10e6, 20e6),
    position_2 = c(20e6, 30e6),
    strand_1 = c("+", "-"),
    strand_2 = c("-", "+"),
    read_count = c(100, 150),
    rearrangement_type = c("deletion", "duplication"),
    type = c("DEL", "DUP")
  )

  # Should NOT give warning
  expect_silent(
    result <- signals:::flip_sv_positions(mock_SV_correct)
  )

  # Data should be unchanged
  expect_equal(result, mock_SV_correct)
})

test_that("flip_sv_positions handles missing columns", {
  # Create mock SV data missing strand columns
  mock_SV_incomplete <- data.frame(
    chromosome_1 = c("11"),
    chromosome_2 = c("11"),
    position_1 = c(20e6),
    position_2 = c(10e6)
    # Missing strand_1, strand_2
  )

  # Should error
  expect_error(
    signals:::flip_sv_positions(mock_SV_incomplete),
    "SV data missing required columns"
  )
})

test_that("plotCNprofile calls flip_sv_positions", {
  # Create mock data with reversed positions
  mock_SV_reversed <- data.frame(
    chromosome_1 = c("11", "11"),
    chromosome_2 = c("11", "11"),
    position_1 = c(20e6, 30e6),
    position_2 = c(10e6, 25e6),
    strand_1 = c("+", "-"),
    strand_2 = c("-", "+"),
    read_count = c(100, 150),
    rearrangement_type = c("deletion", "duplication"),
    type = c("DEL", "DUP")
  )

  # Plotting should trigger the warning
  expect_warning(
    p <- plotCNprofile(mock_CNbins,
                       cellid = "test_cell",
                       SV = mock_SV_reversed,
                       sv_style = "lines_and_arcs",
                       chrfilt = c("11", "6")),
    "Flipped positions and strands for 2 intra-chromosomal SV"
  )

  # Plot should still be created successfully
  expect_true(inherits(p, "ggplot"))
})

test_that("plotCNprofile caps read_count when exceeding sv_read_axis_scale", {
  # Create mock data with high read counts
  mock_SV_high_reads <- data.frame(
    chromosome_1 = c("11", "11"),
    chromosome_2 = c("11", "11"),
    position_1 = c(10e6, 20e6),
    position_2 = c(15e6, 25e6),
    strand_1 = c("+", "-"),
    strand_2 = c("-", "+"),
    read_count = c(500, 1000),  # High read counts
    rearrangement_type = c("deletion", "duplication"),
    type = c("DEL", "DUP")
  )

  # Should give a warning when sv_read_axis_scale is lower than max read_count
  expect_warning(
    p <- plotCNprofile(mock_CNbins,
                       cellid = "test_cell",
                       SV = mock_SV_high_reads,
                       sv_style = "lines_and_arcs",
                       sv_read_axis_scale = 200,
                       chrfilt = c("11", "6")),
    "Capped .* SV breakpoint\\(s\\) with read_count > sv_read_axis_scale"
  )

  # Plot should still be created successfully
  expect_true(inherits(p, "ggplot"))
})

# ---------------------------------------------------------------------------
# SV band above the CN panel (sv_arcs_above = TRUE)
# ---------------------------------------------------------------------------

test_that("classify_sv_side assigns sides by copy number effect", {
  sv <- data.frame(
    chromosome_1 = c("11", "11", "11", "11", "11"),
    chromosome_2 = c("11", "11", "11", "11", "6"),
    position_1   = c(10e6, 20e6, 30e6, 40e6, 50e6),
    position_2   = c(15e6, 25e6, 30e6 + 5e3, 45e6, 10e6),
    strand_1     = c("-", "+", "+", "+", "+"),
    strand_2     = c("+", "-", "+", "+", "-")
  )

  side <- classify_sv_side(sv, rule = "cn_effect", foldback_dist = 30000)
  # -+ duplication = up, +- deletion = down, short ++ foldback = up,
  # long-range ++ inversion = down, translocation = up
  expect_equal(side, c("up", "down", "up", "down", "up"))

  # a foldback is only a foldback within foldback_dist
  expect_equal(classify_sv_side(sv, "cn_effect", foldback_dist = 1000)[3], "down")

  # the other rules
  expect_equal(classify_sv_side(sv, rule = "translocation"),
               c("down", "down", "down", "down", "up"))
  expect_equal(classify_sv_side(sv, rule = "foldback"),
               c("down", "down", "up", "down", "down"))
  expect_error(classify_sv_side(sv, rule = "nonsense"), "sv_arc_side rule must be")
})

test_that("sv_band_trans is monotonic and gives the band a fixed share of the panel", {
  for (base in c("identity", "squashy")) {
    tr <- sv_band_trans(base = base, maxCN = 20, miny = 0,
                        band_frac = 0.3, band_width = 1)
    x <- c(0, 5, 10, 20, 20.5, 21)
    y <- tr$transform(x)

    expect_false(any(is.na(y)))
    expect_true(all(diff(y) > 0))                       # monotonic
    expect_equal(tr$inverse(y), x, tolerance = 1e-6)    # round trips

    # the band is exactly band_frac of the total transformed height
    total <- tr$transform(21) - tr$transform(0)
    band  <- tr$transform(21) - tr$transform(20)
    expect_equal(band / total, 0.3, tolerance = 1e-6)
  }
})

test_that("generate_sv_band_arcs keeps same-bin SVs when min_width is set", {
  # a foldback whose breakends land in one bin has idx_1 == idx_2
  idx1 <- c(10, 30); idx2 <- c(10, 50)

  # without a minimum width the zero-span arc collapses to a single point
  no_min <- generate_sv_band_arcs(idx1, idx2, c("++", "-+"), c("up", "up"),
                                  baseline = 20, half_height = 1, min_width = 0)
  expect_equal(length(unique(no_min$arc_id)), 2)
  expect_equal(diff(range(no_min$idx[no_min$arc_id == "arc_1"])), 0)

  # with one it becomes a narrow but real arc, centred on the breakpoint
  with_min <- generate_sv_band_arcs(idx1, idx2, c("++", "-+"), c("up", "up"),
                                    baseline = 20, half_height = 1, min_width = 4)
  a1 <- with_min[with_min$arc_id == "arc_1", ]
  expect_equal(diff(range(a1$idx)), 4)
  expect_equal(mean(range(a1$idx)), 10)
  expect_gt(max(a1$y), 20)
})

test_that("generate_sv_band_arcs respects side, baseline and the height floor", {
  arcs <- generate_sv_band_arcs(c(10, 10), c(50, 50), c("-+", "+-"),
                                side = c("up", "down"),
                                baseline = 20, half_height = 2, min_frac = 0.15)
  up   <- arcs[arcs$arc_id == "arc_1", ]
  down <- arcs[arcs$arc_id == "arc_2", ]

  expect_true(all(up$y >= 20))     # up arcs stay above the baseline
  expect_true(all(down$y <= 20))   # down arcs stay below it
  # sampled at n_points = 50, which does not land exactly on the t = 0.5 apex
  expect_equal(max(up$y), 22, tolerance = 1e-2)   # longest arc reaches full half height
  expect_equal(min(down$y), 18, tolerance = 1e-2)

  # min_frac sets the floor for the shortest arc
  mixed <- generate_sv_band_arcs(c(10, 10), c(11, 100), c("-+", "-+"),
                                 side = c("up", "up"), baseline = 0,
                                 half_height = 1, min_frac = 0.5)
  short <- max(mixed$y[mixed$arc_id == "arc_1"])
  expect_gte(short, 0.5)
})

test_that("sv_arcs_above reserves the band even when a panel has no SVs", {
  build_y_range <- function(p) {
    b <- ggplot2::ggplot_build(p)
    b$layout$panel_params[[1]]$y.range
  }

  with_sv <- plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                           sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                           sv_arc_side = "cn_effect", chrfilt = c("11", "6"))
  no_sv   <- plotCNprofile(mock_CNbins, cellid = "test_cell",
                           sv_arcs_above = TRUE, sv_arc_side = "cn_effect",
                           chrfilt = c("11", "6"))

  expect_true(inherits(with_sv, "ggplot"))
  expect_true(inherits(no_sv, "ggplot"))
  # both panels must share a y geometry, otherwise stacked panels do not align
  expect_equal(build_y_range(with_sv), build_y_range(no_sv), tolerance = 1e-6)
})

test_that("sv_arcs_above drops the SV read support axis and honours sv_show_lines", {
  p_band <- plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                          sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                          chrfilt = c("11", "6"))
  expect_false(inherits(p_band$scales$get_scales("y")$secondary.axis, "AxisSecondary"))

  # the read count path still builds its secondary axis
  p_reads <- plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                           sv_style = "lines_and_arcs", chrfilt = c("11", "6"))
  expect_true(inherits(p_reads$scales$get_scales("y")$secondary.axis, "AxisSecondary"))

  # turning the lines off removes a layer but keeps the arcs
  p_lines   <- plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                             sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                             sv_show_lines = TRUE, chrfilt = c("11", "6"))
  p_nolines <- plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                             sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                             sv_show_lines = FALSE, chrfilt = c("11", "6"))
  expect_lt(length(p_nolines$layers), length(p_lines$layers))
})

test_that("sv_arcs_above validates its arguments", {
  expect_error(
    plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                  sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                  sv_arc_scale = "nonsense", chrfilt = c("11", "6")),
    "sv_arc_scale must be one of"
  )
  expect_error(
    plotCNprofile(mock_CNbins, cellid = "test_cell", SV = mock_SV,
                  sv_style = "lines_and_arcs", sv_arcs_above = TRUE,
                  sv_band_frac = 1.5, chrfilt = c("11", "6")),
    "sv_band_frac must be a number strictly between 0 and 1"
  )
})

test_that("ybreaks overrides the default y axis breaks", {
  p <- plotCNprofile(mock_CNbins, cellid = "test_cell", chrfilt = c("11", "6"),
                     maxCN = 20, y_axis_trans = "squashy",
                     ybreaks = c(0, 5, 10, 20))
  b <- ggplot2::ggplot_build(p)
  labels <- b$layout$panel_params[[1]]$y$get_labels()
  expect_equal(sort(as.numeric(labels[!is.na(labels)])), c(0, 5, 10, 20))
})
