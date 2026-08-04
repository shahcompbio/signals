library(testthat)
library(ggplot2)
library(dplyr)

# Create mock CNbins data for testing multi-region plotting
set.seed(42)
mock_CNbins_multi_region <- data.frame(
  cell_id = rep("test_cell", 300),
  chr = rep(c("1", "1", "2", "2"), c(100, 100, 50, 50)),
  start = c(
    seq(1, 100e6, by = 1e6),  # Chr1: 1-100 Mb (100 bins)
    seq(150e6, 249e6, by = 1e6),  # Chr1: 150-249 Mb (100 bins)
    seq(1, 50e6, by = 1e6),  # Chr2: 1-50 Mb (50 bins)
    seq(100e6, 149e6, by = 1e6)  # Chr2: 100-149 Mb (50 bins)
  ),
  end = c(
    seq(1e6, 100e6, by = 1e6),
    seq(150e6, 249e6, by = 1e6),
    seq(1e6, 50e6, by = 1e6),
    seq(100e6, 149e6, by = 1e6)
  ),
  state = sample(0:4, 300, replace = TRUE),
  copy = rnorm(300, mean = 2, sd = 0.5)
)

test_that("validate_regions function works correctly", {
  # Valid regions
  valid_regions <- data.frame(
    chr = c("1", "2"),
    start = c(1, 50),
    end = c(10, 75)
  )
  expect_true(validate_regions(valid_regions))
  
  # Invalid: not a data.frame
  expect_error(validate_regions(list(chr = "1", start = 1, end = 10)))
  
  # Invalid: missing columns
  invalid_cols <- data.frame(chr = c("1", "2"), start = c(1, 50))
  expect_error(validate_regions(invalid_cols))
  
  # Invalid: empty data.frame
  expect_error(validate_regions(data.frame()))
  
  # Invalid: start >= end
  invalid_range <- data.frame(
    chr = c("1", "2"),
    start = c(1, 50),
    end = c(10, 40)  # chr2: 50 >= 40
  )
  expect_error(validate_regions(invalid_range))
})

test_that("plotCNprofile with regions parameter works", {
  # Test with single region
  regions_single <- data.frame(
    chr = "1",
    start = 1,
    end = 50
  )
  
  expect_no_error(
    plot <- plotCNprofile(mock_CNbins_multi_region, regions = regions_single)
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile with multiple regions preserves order", {
  # Test with multiple regions in specific order
  regions_ordered <- data.frame(
    chr = c("2", "1", "1"),
    start = c(100, 1, 150),
    end = c(150, 50, 200)
  )
  
  expect_no_error(
    plot <- plotCNprofile(mock_CNbins_multi_region, regions = regions_ordered)
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile with regions and region_gap", {
  regions <- data.frame(
    chr = c("1", "2"),
    start = c(1, 50),
    end = c(30, 100)
  )
  
  # Test with default gap
  expect_no_error(
    plot1 <- plotCNprofile(mock_CNbins_multi_region, regions = regions)
  )
  
  # Test with custom gap
  expect_no_error(
    plot2 <- plotCNprofile(mock_CNbins_multi_region, regions = regions, region_gap = 10)
  )
  
  expect_is(plot1, "ggplot")
  expect_is(plot2, "ggplot")
})

test_that("plotCNprofile with regions and SV data", {
  # Create mock SV data
  mock_SV <- data.frame(
    chromosome_1 = c("1", "1", "2"),
    chromosome_2 = c("1", "2", "2"),
    position_1 = c(10e6, 25e6, 75e6),
    position_2 = c(20e6, 50e6, 100e6),
    strand_1 = c("+", "-", "+"),
    strand_2 = c("-", "+", "-"),
    read_count = c(100, 150, 200),
    rearrangement_type = c("deletion", "translocation", "inversion"),
    type = c("DEL", "TRA", "INV")
  )
  
  regions <- data.frame(
    chr = c("1", "2"),
    start = c(1, 50),
    end = c(30, 100)
  )
  
  # Should work with SV data and curves style
  expect_no_error(
    plot <- plotCNprofile(
      mock_CNbins_multi_region,
      regions = regions,
      SV = mock_SV,
      sv_style = "curves"
    )
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile backward compatibility with chrstart/chrend", {
  # Original functionality should still work
  expect_no_error(
    plot <- plotCNprofile(
      mock_CNbins_multi_region,
      chrfilt = "1",
      chrstart = 1,
      chrend = 50
    )
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile regions takes precedence over chrfilt", {
  regions <- data.frame(
    chr = c("1", "2"),
    start = c(1, 50),
    end = c(30, 100)
  )
  
  # regions should override chrfilt
  expect_no_error(
    plot <- plotCNprofile(
      mock_CNbins_multi_region,
      chrfilt = "1",  # This should be ignored
      regions = regions
    )
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile handles same chromosome regions", {
  # Multiple regions from same chromosome
  regions <- data.frame(
    chr = c("1", "1"),
    start = c(1, 150),
    end = c(50, 200)
  )
  
  expect_no_error(
    plot <- plotCNprofile(
      mock_CNbins_multi_region,
      regions = regions,
      region_gap = 5
    )
  )
  expect_is(plot, "ggplot")
})

test_that("plotCNprofile with returnlist = TRUE works with regions", {
  regions <- data.frame(
    chr = "1",
    start = 1,
    end = 50
  )
  
  result <- plotCNprofile(
    mock_CNbins_multi_region,
    regions = regions,
    returnlist = TRUE
  )
  
  expect_is(result, "list")
  expect_named(result, c("CN", "plist"))
  expect_is(result$CN, "ggplot")
  expect_is(result$plist, "list")
})

test_that("plotCNprofile with no data in regions returns error", {
  # Regions outside of the data range
  regions <- data.frame(
    chr = "10",  # This chromosome doesn't exist in the data
    start = 1,
    end = 50
  )
  
  expect_error(
    plotCNprofile(mock_CNbins_multi_region, regions = regions)
  )
})
