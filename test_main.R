#!/usr/bin/Rscript
source("main.R")
library(testthat)

# first argument is just a description, second is a code block inside {}
test_that("read_data loads correctly", {
  # since we're sourcing main.R, we can just call the function from there as is.
  # in this case we can use the intensity data, but if this file is changed the 
  # test may break. Store results into res.
  res <- read_data("example_intensity_data.csv", " ")
  # first test is checking the dimensions. expect equal passes if the two 
  # arguments are the same, so we can just make sure it has the right number of 
  # rows and columns. Also good for students to see what's expected (hopefully).
  expect_equal(dim(res), c(54675, 35))
  # since you specify a data.frame and not a tibble in the function description, 
  # we can check the class() to make sure it is a df and not anything else. 
  expect_equal("data.frame", class(res))
  # this is probably enough tests for a simple data loading function like this
})

test_that("calculate_variance", {
  skip(message = "Implement later.")
})

test_that("make_variance_tibble", {
  skip(message = "Implement later.")
})

test_that("plot_pca_variance creates a ggplot object with the correct geoms", {
  pcs <- factor(paste0("PC", 1:35), levels = paste0("PC", 1:35), ordered = T)
  tibble(variance_explained = runif(35, 1.0e-30, 1.5e-01),
         principal_components = pcs) %>%
    mutate(cumulative = cumsum(variance_explained)) -> variance_tibble
  plot <- plot_pca_variance(variance_tibble)
  
  # testing all geoms in the ggplot object
  geoms <- c()
  for (geom in plot$layers) {
    geoms <- c(geoms, class(geom$geom)[1])
  }
  expect_true("GeomBar" %in% geoms)
  expect_true("GeomLine" %in% geoms)
  expect_true("GeomPoint" %in% geoms)
})

test_that("make_biplot", {
  skip(message = "Implement later.")
})

test_that("list_significant_probes", {
  skip(message = "Implement later.")
})

test_that("return_de_intensity", {
  skip(message = "Implement later.")
})

test_that("plot_heatmap has correct rows and cols", {
  # this one is an interesting case, since its object (what the function returns)
  # doesn't have as many descriptive elements as ggplot's. We can at least 
  # confirm that their function creates a heatmap and that it has a predictable 
  # number of rows/columns. Since we don't know if their earlier functions will 
  # work, we will create a fake de_intensity matrix to test with.
  fx_intensity <- data.frame(GSM1 = runif(1000, 0, 15),
                             GSM2 = runif(1000, 0, 15),
                             GSM3 = runif(1000, 0, 15),
                             GSM4 = runif(1000, 0, 15))
  row.names(fx_intensity) <- paste0(1:1000, "_at")
  fx_intensity <- as.matrix(fx_intensity)
  heatmap <- plot_heatmap(fx_intensity)
  expect_equal(length(heatmap$colInd), 4)
  expect_equal(length(heatmap$rowInd), 1000)
  expect_equal(class(heatmap), "list")
})
