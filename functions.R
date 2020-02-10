########## Functions ##########

# Odds ratio
odds.ratio <- function(a, b, total){
  a.b <- length(intersect(a, b))
  a <- length(a) 
  b <- length(b) 
  a.nonb <- a - a.b
  nona.b <- b - a.b
  nona.nonb <- total - a.b - a.nonb - nona.b
  v <- c(a.b, a.nonb, nona.b, nona.nonb)
  x <- matrix(v, 2, byrow = TRUE)
  or <- OddsRatio(x) # Odds of b
  or
}

# Hypergeometric test
hyper.test <- function(a, b, total){
  genes <- intersect(a, b)
  overlap <- length(genes)
  ns1 <- length(a)
  ns2 <- length(b)
  p <- phyper(overlap - 1, ns1, total - ns1, ns2, lower.tail = FALSE)
  p
}

# Test for two lists of gene sets and correct P for cell-types tested
hypertest.or.table <- function(l1, l2, total){ # two lists of gene sets
  t <- lapply(l1, function(set1){
    # Overlap with each gene set
    t <- t(sapply(l2, function(set2){
      or <- odds.ratio(set2, set1, total) # odds of l2
      p <- hyper.test(set2, set1, total)
      c(or = or, pvalue = p)
    }))
    bh <- p.adjust(t[,"pvalue"], method = "BH")
    cbind(t, bh) # 2D-table: l2 variables x measures
  })
  simplify2array(t) # 3d-array: l2 variables x measures x l1 variables
}