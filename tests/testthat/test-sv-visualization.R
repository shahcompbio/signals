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

  # Dimensions should be different for horizontal vs vertical
  expect_false(identical(dim(leg_vert), dim(leg_horiz)))
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
