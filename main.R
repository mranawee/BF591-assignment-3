read_data <- function(intensity_data, delimiter) {
  data <- as.data.frame(read.csv(intensity_data, sep= delimiter))
  return(data)
}

#Calculate PCA variance explained
pca_var <- function(pca_results) {
  ve <- pca_results$sdev^2
  pve <- ve / sum(ve)
  return(pve)
}

make_var_tib <- function(pca_var, pca_results) {
  var_tib <- tibble(variance_explained = pca_var) %>% 
  mutate(principal_components = factor(colnames(pca_results$x), levels=colnames(pca_results$x)), cumulative = cumsum(variance_explained))
  return(var_tib)
}

plot_pca_ve <- function(var_tib) {
  ls_cols <- c('Variance Explained'='black')
  bar_cols <-c('Cumulative'='light blue')
  var_tib %>%
    ggplot() + 
    geom_bar(aes(x=principal_components, y=variance_explained, fill='Variance Explained'), stat='identity', color='black') + 
    geom_line(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) +
    geom_point(aes(x=principal_components, y=cumulative, group=1, color='Cumulative')) +
    scale_colour_manual(name="Cumulative",values=ls_cols) +
    scale_fill_manual(name="Variance Explained",values=bar_cols) +
    labs(x='PC', y='% variance') +
    theme_classic(base_size=8) + 
    theme(axis.text.x=element_text(angle=90,hjust=1))
  
}

add_pca_labels <- function(metadata, pca_results) {
  meta <- readr::read_csv(metadata)
  class <- meta %>% select(geo_accession, SixSubtypesClassification)
  x <- pca_results$x %>% as_tibble(rownames='geo_accession') %>% left_join(class, by='geo_accession')
  return(x)
}

plot_pca <- function(labeled) {
  biplot <- labeled %>% ggplot() + geom_point(aes(x=PC1, y=PC2, color=SixSubtypesClassification))
  return(biplot)
}

list_sig_ids <- function(diff_exp_csv, fdr_threshold) {
sig_ids <- read_data('differential_expression_results.csv', ',') %>% 
  as_tibble(rownames='probeid') %>% 
  filter(padj < fdr_threshold) %>% 
  select(probeid) %>% 
  as.list()
  return(sig_ids$probeid)
}

return_de_intensity <- function(intensity, sig_ids_list) {
  de_intensity <- intensity %>% 
    as_tibble(rownames='probeid') %>% 
    filter(probeid %in% sig_ids_list) %>% 
    column_to_rownames(var='probeid') %>% 
    as.matrix()
  return(de_intensity)
}

plot_heatmap <- function(scaled, num_cols, color_pal) {
  col.pal <- RColorBrewer::brewer.pal(num_cols, color_pal)
  return(heatmap(de_intensity, col=col.pal))
}