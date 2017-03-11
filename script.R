library(dplyr)
library(limma)

res.data <- read.csv("data.csv")

fit <- factanal(res.data, factors=8, rotation="promax", trace=T)
# unable to optimize from this starting value
fit <- factanal(res.data, factors=7, rotation="promax", trace=T)
# unable to optimize from this starting value

## Features which make factanal fail

factanal_fail_feats <- function(data, proportion=.75, n=10000, mod=1000) {
  col_tot <- ncol(data)
  col_sel <- col_tot * proportion
  feat_scores <- rep(0, col_tot)
  names(feat_scores) <- colnames(data)
  for (i in 1:n) {
    selected <- sample(1:col_tot, col_sel)
    idata <- data[, selected]
    tryCatch({
      factanal(idata, factors=8, rotation="promax")
    },
    error=function(e) {
      feat_scores[selected] <<- feat_scores[selected] + 1
    })
    if (i %% mod == 0) {
      print(head(sort(feat_scores, decreasing=TRUE)))
    }
  }
  sort(feat_scores, decreasing=TRUE)
}

# WARNING: Takes a while to complete! Interim results are reported periodically.
fail_feats <- factanal_fail_feats(res.data)
head(fail_feats)
#   NN  NOM   VF STA4  GEN CLUN
# 4459 4406 4190 4186 4088 4063
fit <- factanal(select(res.data, -NN), factors=8, rotation="promax", trace=T)
# unable to optimize from this starting value
fit <- factanal(select(res.data, -NN, -NOM), factors=8, rotation="promax", trace=T)
# unable to optimize from this starting value
fit <- factanal(select(res.data, -NN), factors=7, rotation="promax", trace=T)
# success!

## Features which cause linear dependence between columns

lin_dep_feats <- function(data, proportion=.75, n=10000) {
  col_tot <- ncol(data)
  col_sel <- col_tot * proportion
  feat_scores <- rep(0, col_tot)
  names(feat_scores) <- colnames(data)
  for (i in 1:n) {
    selected <- sample(1:col_tot, col_sel)
    idata <- data[, selected]
    # is.fullrank checks if matrix has full column rank; nonEstimable returns column identifiers of
    # columns which are linearly dependent on previous ones; both come from the limma package
    if (!is.fullrank(idata)) feat_scores[selected] <- feat_scores[selected] + 1
  }
  sort(feat_scores, decreasing=TRUE)
}

is.fullrank(res.data)
# FALSE
ld_feats <- lin_dep_feats(res.data)
head(ld_feats)
#   WL   TC NEG1 COH1  NNV VTE1
# 7505 7386 7065 7062 7057 7056
is.fullrank(select(res.data, -WL, -TC))
# TRUE
fit <- factanal(select(res.data, -WL, -TC), factors=8, rotation="promax", trace=T)
# unable to optimize from this starting value
fit <- factanal(select(res.data, -WL, -TC), factors=7, rotation="promax", trace=T)
# unable to optimize from this starting value

# If I remove the 4 features identified as "problematic" by both approaches, the factanal call
# succeeds even with 8 factors.
fit <- factanal(select(res.data, -NN, -NOM, -WL, -TC), factors=8, rotation="promax", trace=T)
# success!

# Randomly removing 4 features doesn't do anything helpful (whew).
fit <- factanal(res.data[, sample(1:ncol(res.data), ncol(res.data)-4)], factors=8, rotation="promax", trace=T)
# unable to optimize from this starting value

## Correlations between "problematic" features and the rest

feat_cor <- function(data, feat_name) {
  cor <- cor(data[[feat_name]], data)
  sort(setNames(c(abs(cor)), colnames(cor)), decreasing=TRUE)
}

head(feat_cor(res.data, "NN"))
#        NN      CLUN        PP       GEN        WL        DB
# 1.0000000 0.9172335 0.9017139 0.8935881 0.8575989 0.8565588
head(feat_cor(res.data, "WL"))
#        WL      ATA1        AA       GEN      STA4        NN
# 1.0000000 0.9053217 0.8816113 0.8799450 0.8589837 0.8575989
head(feat_cor(res.data, "TC"))
#        TC       PPI       BIG       NNV       VVO      STA4
# 1.0000000 0.5262275 0.4937046 0.4823765 0.4687399 0.4630486

# Nothing too extreme here, right?
