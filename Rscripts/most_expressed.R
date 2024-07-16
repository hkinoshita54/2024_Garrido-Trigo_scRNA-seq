C <- seu[["RNA"]]$counts
# C@x <- C@x / rep.int(colSums(C), diff(C@p)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])),
        cex = 0.1, las = 1, xlab = "Percent counts per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE
)
