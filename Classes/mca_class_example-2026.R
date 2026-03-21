
# Analysis ===========================
# Inertia in indicator
survey_data <- read.table("data/survey_example.dat")
survey_data[] <- lapply(survey_data, as.factor)
indicator_matrix <- tab.disjonctif(survey_data)
print(head(indicator_matrix))

burt_matrix <- t(indicator_matrix) %*% indicator_matrix
print(head(burt_matrix))

Q <- ncol(survey_data)
J <- sum(sapply(survey_data, nlevels))

(I_I <- J/Q -1)
# Inertia in Burt
n_k <- diag(burt_matrix)

expected_denominators <- outer(n_k, n_k)

(I_B_direct <- (1 / Q**2) * sum((burt_matrix**2) / expected_denominators) - 1)
# Adjusted Burt:
(I_adj <- 1/(Q**2-Q)*(Q**2*I_B_direct-(J-Q)))

# MCA on Indicator
res.mca.ind <- MCA(survey_data, graph=FALSE, method="Indicator")
res.mca.ind$eig
round(res.mca.ind$eig[,1], 3)
round(res.mca.ind$svd$vs, 3)
round(res.mca.ind$svd$vs, 3)**2
cumsum(res.mca.ind$eig[,1])
# adjusted inertias
filter.eig <- (res.mca.ind$eig[,1])> 1/Q
filtered.eigen <- res.mca.ind$eig[,1][filter.eig]
adj.eig <- ((Q/(Q-1))**2)*((filtered.eigen)-1/Q)**2
sum(adj.eig)
cumsum(adj.eig)
df.adj.eig <- data.frame(
  values = adj.eig,
  percent = adj.eig/sum(adj.eig)*100,
  cumul = cumsum(adj.eig)
)

# MCA on Burt
res.mca.burt <- MCA(survey_data, graph=FALSE, method="Burt")
res.mca.burt$eig
round(res.mca.burt$eig[,1], 3)
round(res.mca.burt$svd$vs, 3)
round(res.mca.burt$svd$vs, 3)**2
cumsum(res.mca.burt$eig[,1])

# adjusted inertias
filter.eig <- sqrt(res.mca.burt$eig[,1])> 1/Q
filtered.eigen <- res.mca.burt$eig[,1][filter.eig]
adj.eig <- ((Q/(Q-1))**2)*(sqrt(filtered.eigen)-1/Q)**2
sum(adj.eig)
cumsum(adj.eig)
df.adj.eig <- data.frame(
  values = adj.eig,
  percent = adj.eig/sum(adj.eig)*100,
  cumul = cumsum(adj.eig)
)

# Biplots
## FactoMineR by default produces a symmetric biplot.
## We can interpret distances between individuals and distances between variables
## BUT NEVER distances between individuals and variables. BUT, the direction in
## the space provides interpretation about individuals WITH the variables. This
## is the BARYCENTRIC RULE and it means that categories pull individuals in their
## direction.
g1 <- fviz_mca(res.mca.ind, title="MCA-Biplot Indicator")
g2 <- fviz_mca(res.mca.burt, title="MCA-Biplot Burt")
g1
g2

fviz_mca_biplot(res.mca.ind, col.var = "cos2", select.ind= list(cos2=25), repel = TRUE)
fviz_mca_biplot(res.mca.ind, col.ind = "cos2", select.ind= list(cos2=25), repel = TRUE)

## Interpretations ============
### eta2 relates questions with dimensions (factor axes)
### v.test relates direction (sign) with answers when |magnitude|>1.96
summary(res.mca.ind)
summary(res.mca.burt)
