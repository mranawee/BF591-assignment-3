---
title: "Untitled"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
test_sd <- sort(runif(20, min=0, max=100), decreasing=TRUE)
test_sd = list(sdev = test_sd)
sd <- test_sd$sdev
```

```{r}
calculate_variance_explained <- function(pca_results) {
  v <- pca_results$sdev^2
  ve <- v / sum(v)
  return(ve)
}


```



```{r}




#v <- sd^2
#ve <- v / sum(v)

test_that('function calculates variance and ve', {

  test_variance <- test_sd$sdev^2
  test_ve <- test_variance / sum(test_variance)
  
  expect_equal(calculate_variance_explained(test_sd), test_ve)
})
  
  


```


```{r}

```