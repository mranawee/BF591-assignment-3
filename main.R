library('tidyverse')
library('RColorBrewer')
library('readr')
library(ggplot2)
library(grid)
library(gridExtra)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("sinhrks/ggfortify")
library(ggfortify)

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
  df <- read.delim(intensity_data, header = TRUE,sep = delimiter )
  return(df)
}

intensity <- read_data('example_intensity_data.csv', ' ')
head(intensity, 2)
dim(intensity)
#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
#' 

pca_results <- prcomp(scale(t(intensity)), center=FALSE, scale=FALSE)

calculate_variance_explained <- function(pca_results){
  return(pca_results$sdev^2 / sum(pca_results$sdev^2))
}


#pca_results$var <- (pca_results$sdev)^2
#pca_results$var

pc_variance_explained <- calculate_variance_explained(pca_results)
pc_variance_explained

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  var_tibble<- as_tibble(pca_ve)
  colnames(var_tibble) <- 'variance_explained'
  var_tibble <- mutate(var_tibble, principal_components = colnames(pca_results$rotation))
  var_tibble$principal_components <- as.factor(var_tibble$principal_components)
  var_tibble <- mutate(var_tibble, cumulative = cumsum(variance_explained))
  
  return(var_tibble)
}

variance_tibble <- make_variance_tibble(pc_variance_explained, pca_results)
head(variance_tibble, 5)

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble gnerated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples

plot_pca_variance <- function(variance_tibble) {
  plt <- ggplot(variance_tibble) +
    geom_bar(stat = 'identity', aes(x = reorder(principal_components, -variance_explained),
                                    y = variance_explained,
                                    fill = "Variance Explained")) +
    geom_point(aes(x=principal_components, y=cumulative, group = 1, colour = "Cumulative")) +
    geom_line(aes(x=principal_components, y=cumulative, group = 1, colour = "Cumulative")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_colour_manual(name="Cumulative",values='black') +
    scale_fill_manual(name="Variance Explained",values='light blue') 
    return(plt + xlab('PC') +ylab('% variance'))
}

plot_pca_variance(variance_tibble)


  #geom_line(aes(x= principal_components, y=cumulative))

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples

make_biplot <- function(metadata, pca_results) {
  meta <- read_csv(metadata) 
  df <- as.data.frame(pca_results$x)
  df$SixSubtypesClassification <- meta$SixSubtypesClassification[100:134]
  p <- ggplot(df, aes(x=PC1, y=PC2, color=SixSubtypesClassification)) + geom_point()
  return(p)
}

make_biplot('proj_metadata.csv', pca_results)


#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_csv (str): The path to the differential expression results
#'   file we have provided
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the CSV.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
  data <- read_data(diff_exp_csv,',')
  filtered <- filter(data, padj < fdr_threshold)
  return(rownames(filtered))
}

sig_ids <- list_significant_probes('differential_expression_results.csv', .01)
head(sig_ids)
length(sig_ids)

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  de_intensity <- filter(intensity, rownames(intensity) %in% sig_ids_list)
  return(as.matrix(de_intensity))
}

de_intensity <- return_de_intensity(intensity, sig_ids)
typeof(head(de_intensity, 2))

head(intensity)
de_intensity <- filter(intensity, rownames(intensity) %in% sig_ids)
head(de_intensity,2)

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  col.pal <- brewer.pal(num_colors, palette)
  map <- heatmap(de_intensity, col= col.pal)
  return(map)
}

plot_heatmap(de_intensity, 11, 'RdBu')

#ggplot 4
#barplot is geom
#line is geom point, geom bar, geom line

#for 5, use geom point
