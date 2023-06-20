############################################################################
#                       Initialize working space                           #
############################################################################
# 1. remove all objects to clean the working space
rm(list = ls())
# 2. show the current working directory
print(sprintf("current working directory: %s", getwd()))
# 3. set random seed
set.seed(22200056)

############################################################################
#                                Question 01                               #
############################################################################
# 1. load dataset, (431, 582), global_data is shared among all questions
global_data <- read.csv("dataset/Milk_MIR_Traits_data_2023.csv")
# 2. remove a random observation, (430, 582)
q1.to_remove <- sample(1:nrow(global_data), 1)
global_data <- global_data[-q1.to_remove, ]

############################################################################
#                                Question 02                               #
############################################################################
# 1. remove observations with NA beta_lactoglobulin_b values, (305, 582)
global_data <- global_data[!is.na(global_data$beta_lactoglobulin_b), ]

# 2. spectra visualization
# extract spectra data (the last 531 columns)
q2.spectra <- as.matrix(global_data[, c(52:582)])
# extract column names (wavelengths) as labels of x-axis
q2.wavelen <- colnames(q2.spectra)
for (i in 1:531) {
    q2.wavelen[i] <- substr(q2.wavelen[i], 2, nchar(q2.wavelen[i]))
}
# line graph for spectra data
# matplot: plot columns of matrix
matplot(
    t(q2.spectra),
    type = "l",
    xlab = "wavelength (nm)",
    ylab = "absorbance",
    axes = FALSE
)
q2.x_inds <- seq(1, ncol(q2.spectra), length = 20)
axis(1, at = q2.x_inds, labels = q2.wavelen[q2.x_inds])
axis(2)
# 3. beta_lactoglobulin_b visualization
# extract beta_lactoglobulin_b data
q2.beta_ln_b <- global_data$beta_lactoglobulin_b
# histogram for beta_lactoglobulin_b data
hist(
    q2.beta_ln_b,
    xlab = "beta_lactoglobulin_b (gram)",
    ylab = "frequency",
    breaks = 100,
    col = "lightblue"
)
# pie chart for beta_lactoglobulin_b
q2.piechart.x <- c(sum(q2.beta_ln_b == 0.0), sum(q2.beta_ln_b != 0.0))
q2.piechart.percent <- round(100 * q2.piechart.x / sum(q2.piechart.x), 1)
pie(
    q2.piechart.x,
    labels = q2.piechart.percent,
    col = rainbow(length(q2.piechart.x))
)
legend(
    "topright",
    c("zero", "non-zero"),
    cex = 0.8,
    fill = rainbow(length(q2.piechart.x))
)

# 3. handle outliers
# remove observations with beta_lactoglobulin_b outside of 3 standard
# deviations from mean, (303, 582)
# upper = mean + 3 * std
q2.upper_bound <- mean(q2.beta_ln_b) + 3 * sd(q2.beta_ln_b)
# lower = mean - 3 * std
q2.lower_bound <- mean(q2.beta_ln_b) - 3 * sd(q2.beta_ln_b)
q2.to_remove <-
    which(q2.beta_ln_b < q2.lower_bound | q2.beta_ln_b > q2.upper_bound)
global_data <- global_data[-q2.to_remove, ]

############################################################################
#                                Question 03                               #
############################################################################
q3.spectra <- as.matrix(global_data[, c(52:582)])

# calculate dissimilarity matrix, Euclidean distance is employed
q3.dist_mat <- dist(q3.spectra, method = "euclidean")

# 1. hierarchical clustering
# single linkage
q3.hier_single <- hclust(q3.dist_mat, method = "single")
plot(q3.hier_single)
# complete linkage
q3.hier_complete <- hclust(q3.dist_mat, method = "complete")
plot(q3.hier_complete)
# average linkage
q3.hier_average <- hclust(q3.dist_mat, method = "average")
plot(q3.hier_average)

# choose complete linkage and partition data into 3 and 5 clusters
q3.hier_complete_3 <- cutree(q3.hier_complete, 3)
print(table(q3.hier_complete_3))

q3.hier_complete_5 <- cutree(q3.hier_complete, 5)
print(table(q3.hier_complete_5))

# 2. k-means clustering
q3.ss <- rep(0, 10)
# calculate within group sum of squares for k = 1
# (explanation is given in the report)
q3.ss[1] <- (nrow(global_data) - 1) * sum(apply(q3.spectra, 2, var))
# calculate within group sum of squares for k = 2...10
for (k in 2:10) {
    q3.ss[k] <- sum(kmeans(q3.spectra, centers = k)$withinss)
}
plot(1:10, q3.ss, type = "b", xlab = "k", ylab = "SS")

# partition data into 3 and 5 clusters
q3.kmeans_3 <- kmeans(q3.spectra, centers = 3)
print(table(q3.kmeans_3$cluster))
q3.kmeans_5 <- kmeans(q3.spectra, centers = 5)
print(table(q3.kmeans_5$cluster))

# 3. compare hierarchical clustering with k-means clustering
library(e1071)
print(table(q3.hier_complete_3, q3.kmeans_3$cluster))
print(classAgreement(table(q3.hier_complete_3, q3.kmeans_3$cluster))$crand)
print(table(q3.hier_complete_5, q3.kmeans_5$cluster))
print(classAgreement(table(q3.hier_complete_5, q3.kmeans_5$cluster))$crand)

# 4. compare clustering solution with the known origins
# Breed
print(table(q3.hier_complete_3, global_data$Breed))
print(classAgreement(table(q3.hier_complete_3, global_data$Breed))$crand)

# Parity
print(table(q3.hier_complete_3, global_data$Parity))
print(classAgreement(table(q3.hier_complete_3, global_data$Parity))$crand)

# Milking_Time
print(table(q3.hier_complete_3, global_data$Milking_Time))
print(classAgreement(table(q3.hier_complete_3, global_data$Milking_Time))$crand)

# DaysInMilk
q3.upper_quantile <- quantile(global_data$DaysInMilk, probs = 2 / 3)
q3.lower_quantile <- quantile(global_data$DaysInMilk, probs = 1 / 3)
q3.is_high <- global_data$DaysInMilk >= q3.upper_quantile
q3.is_medium <- global_data$DaysInMilk >= q3.lower_quantile &
    global_data$DaysInMilk < q3.upper_quantile
q3.is_low <- global_data$DaysInMilk < q3.lower_quantile
q3.days_in_milk <- rep(NA, nrow(global_data))
q3.days_in_milk[q3.is_high] <- 3
q3.days_in_milk[q3.is_medium] <- 2
q3.days_in_milk[q3.is_low] <- 1
print(table(q3.hier_complete_3, q3.days_in_milk))
print(classAgreement(table(q3.hier_complete_3, q3.days_in_milk))$crand)

############################################################################
#                                Question 04                               #
############################################################################
q4.spectra <- as.matrix(global_data[, c(52:582)])

# apply PCA to spectra data
q4.pca <- prcomp(q4.spectra)
# variance of each PC
q4.var <- q4.pca$sdev^2
# proportion of variance explained
q4.prop <- q4.var / sum(q4.var)
# cumulative proportion of variance explained
q4.cum_prop <- cumsum(q4.prop)
print(q4.cum_prop[1:10])
# plot cumulative proportion of variance explained
plot(q4.cum_prop[1:10],
    xlab = "principal component",
    ylab = "cumulative proportion of variance explained",
    ylim = c(0.6, 1.0),
    type = "b"
)

############################################################################
#                                Question 05                               #
############################################################################
q5.spectra <- as.matrix(global_data[, c(52:582)])

# apply PCA to spectra data
q5.pca <- prcomp(q5.spectra)
# only use loadings for the first 3 PCs
q5.loadings_3 <- q5.pca$rotation[, 1:3]
# center spectra data, (303, 531)
q5.centered_spectra <- scale(
    q5.spectra,
    center = TRUE,
    scale = FALSE
)
# scores on the first 3 PCs for spectra data, (303, 3)
q5.scores_3 <- q5.centered_spectra %*% q5.loadings_3

# compare results with
# 1) builtin scores from result of prcomp command
# 2) scores from predict command
# to check the correctness
q5.scores_3_prcomp <- q5.pca$x[, 1:3]
stopifnot(q5.scores_3 == q5.scores_3_prcomp)
q5.scores_3_predict <- predict(q5.pca)[, 1:3]
stopifnot(q5.scores_3 == q5.scores_3_predict)

# scatter plot for all pairs of PCs
pairs(q5.scores_3)

# differentiate the observations according to hierarchical clustering
pairs(q5.scores_3, col = q3.hier_complete_3, pch = q3.hier_complete_3)

q5.colors <- c(
    "black",
    "blue",
    "orange",
    "yellow",
    "red",
    "green",
    "grey",
    "tomato",
    "skyblue",
    "cyan"
)
# differentiate the observations according to Breed
q5.breed.cols <- rep(NA, nrow(global_data))
q5.breeds <- unique(global_data$Breed)
for (i in 1:length(q5.breeds)) {
    q5.breed.cols[global_data$Breed == q5.breeds[i]] <- q5.colors[i]
}
pairs(q5.scores_3, col = q5.breed.cols)

# differentiate the observations according to Parity
q5.parity.cols <- rep(NA, nrow(global_data))
q5.parities <- unique(global_data$Parity)
for (i in 1:length(q5.parities)) {
    q5.parity.cols[global_data$Parity == q5.parities[i]] <- q5.colors[i]
}
pairs(q5.scores_3, col = q5.parity.cols)

# differentiate the observations according to Milking_Time
q5.milking_time.cols <- rep(NA, nrow(global_data))
q5.milking_times <- unique(global_data$Milking_Time)
for (i in 1:length(q5.milking_times)) {
    q5.milking_time.cols[global_data$Milking_Time == q5.milking_times[i]] <- q5.colors[i]
}
pairs(q5.scores_3, col = q5.milking_time.cols)

# differentiate the observations according to DaysInMilk
q5.upper_quantile <- quantile(global_data$DaysInMilk, probs = 2 / 3)
q5.lower_quantile <- quantile(global_data$DaysInMilk, probs = 1 / 3)
q5.is_high <- global_data$DaysInMilk >= q5.upper_quantile
q5.is_medium <- global_data$DaysInMilk >= q5.lower_quantile &
    global_data$DaysInMilk < q5.upper_quantile
q5.is_low <- global_data$DaysInMilk < q5.lower_quantile
q5.days_in_milk <- rep(NA, nrow(global_data))
q5.days_in_milk[q5.is_high] <- 3
q5.days_in_milk[q5.is_medium] <- 2
q5.days_in_milk[q5.is_low] <- 1
q5.days_in_milk.cols <- rep(NA, nrow(global_data))
q5.days_in_milks <- unique(q5.days_in_milk)
for (i in 1:length(q5.days_in_milks)) {
    q5.days_in_milk.cols[q5.days_in_milk == q5.days_in_milks[i]] <- q5.colors[i]
}
pairs(q5.scores_3, col = q5.days_in_milk.cols)

############################################################################
#                                Question 07                               #
############################################################################
library(pls)

# split dataset into train set and test set
q7.n_observations <- nrow(global_data)
# train set
q7.train_inds <- sample(
    1:q7.n_observations,
    size = 2 / 3 * q7.n_observations,
    replace = FALSE
)
q7.train <- data.frame(
    spectra = I(as.matrix(global_data[q7.train_inds, c(52:582)])),
    beta_ln_b = global_data[q7.train_inds, "beta_lactoglobulin_b"]
)
# test set
q7.test_inds <- c(1:q7.n_observations)[-q7.train_inds]
q7.test <- data.frame(
    spectra = I(as.matrix(global_data[q7.test_inds, c(52:582)])),
    beta_ln_b = global_data[q7.test_inds, "beta_lactoglobulin_b"]
)

# PCR fitting
q7.pcr <- pcr(
    formula = q7.train$beta_ln_b ~ q7.train$spectra,
    ncomp = 100,
    data = q7.train,
    center = TRUE,
    scale = FALSE,
    validation = "CV"
)
print(summary(q7.pcr))
# plot train observations and regression line
plot(
    RMSEP(q7.pcr),
    legendpos = "topleft",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
plot(q7.pcr, ncomp = 37, asp = 1, line = TRUE)
# PCR prediction
q7.rmse <- rep(NA, 100)
for (i in 1:100) {
    q7.preds <- predict(q7.pcr, ncomp = i, newdata = q7.test$spectra)
    q7.rmse[i] <- sqrt(mean((q7.preds - q7.test$beta_ln_b)^2))
}
plot(
    q7.rmse,
    type = "l",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
print(min(q7.rmse))

############################################################################
#                                Question 08                               #
############################################################################
impute <- function(X, method = "pca", k = 1, tol = 1e-7, maxiter = 1000) {
    # step 1: fill missing values with corresponding column mean
    Xtilde <- X
    xbar <- colMeans(X, na.rm = TRUE)
    for (col_idx in seq_len(ncol(Xtilde))) {
        na_inds <- which(is.na(Xtilde[, col_idx]))
        Xtilde[na_inds, col_idx] <- xbar[col_idx]
    }
    if (method == "mean") {
        return(Xtilde)
    }
    # step 2: repeat step 2(a) to step 2(c) until the objective doesn't decrease
    isna <- is.na(X)
    mse_prev <- mean((scale(X, center = xbar, scale = FALSE)[!isna])^2)
    for (iter in 1:maxiter) {
        # step 2(a): recoustruct Xtilde
        if (method == "pca") {
            pca <- prcomp(Xtilde)
            scores_k <- pca$x[, 1:k]
            loadings_k <- pca$rotation[, 1:k]
            Xhat <- scores_k %*% t(loadings_k)
            Xhat <- scale(
                Xhat,
                center = -colMeans(Xtilde),
                scale = FALSE
            )
        } else if (method == "svd") {
            svd_ret <- svd(Xtilde)
            Xhat <- svd_ret$u[, 1:k] %*% (svd_ret$d[1:k] * t(svd_ret$v[, 1:k]))
        } else {
            stopifnot(FALSE)
        }
        # Step 2(b): fill missing values with corresponding values in Xhat
        Xtilde[isna] <- Xhat[isna]
        # Step 2(c): calculate mse and relative error
        # calculate mse
        mse <- mean(((X - Xhat)[!isna])^2)
        # calculate relative error
        rel_err <- (mse_prev - mse) / mean(X[!isna]^2)
        # update mse_prev
        mse_prev <- mse
        cat("Iter:", iter, "mse:", mse, "Rel. Err:", rel_err, "\n")
        # termination condition
        if (rel_err <= tol) {
            break
        }
    }
    return(Xtilde)
}

# test impute function on lena128 dataset
library(filling)
data(lena128)
q8.lena128 <- aux.rndmissing(lena128, x=0.05)

# mean
fill1_mean <- as.matrix(impute(q8.lena128, method="mean", k=10, tol = 1e-3))
fill2_mean <- as.matrix(impute(q8.lena128, method="mean", k=25, tol = 1e-3))
fill3_mean <- as.matrix(impute(q8.lena128, method="mean", k=50, tol = 1e-3))
opar_mean <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
image(q8.lena128, col=gray((0:100)/100), axes=FALSE, main="5% missing")
image(fill1_mean, col=gray((0:100)/100), axes=FALSE, main="10 regressors")
image(fill2_mean, col=gray((0:100)/100), axes=FALSE, main="25 regressors")
image(fill3_mean, col=gray((0:100)/100), axes=FALSE, main="50 regressors")
par(opar_mean)
mtext("mean", side = 3, line = -21, outer = TRUE)

# pca
fill1_pca <- as.matrix(impute(q8.lena128, method="pca", k=10, tol = 1e-3))
fill2_pca <- as.matrix(impute(q8.lena128, method="pca", k=25, tol = 1e-3))
fill3_pca <- as.matrix(impute(q8.lena128, method="pca", k=50, tol = 1e-3))
opar_pca <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
image(q8.lena128, col=gray((0:100)/100), axes=FALSE, main="5% missing")
image(fill1_pca, col=gray((0:100)/100), axes=FALSE, main="10 regressors")
image(fill2_pca, col=gray((0:100)/100), axes=FALSE, main="25 regressors")
image(fill3_pca, col=gray((0:100)/100), axes=FALSE, main="50 regressors")
par(opar_pca)
mtext("pca", side = 3, line = -21, outer = TRUE)

# svd
fill1_svd <- as.matrix(impute(q8.lena128, method="svd", k=10, tol = 1e-3))
fill2_svd <- as.matrix(impute(q8.lena128, method="svd", k=25, tol = 1e-3))
fill3_svd <- as.matrix(impute(q8.lena128, method="svd", k=50, tol = 1e-3))
opar_svd <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
image(q8.lena128, col=gray((0:100)/100), axes=FALSE, main="5% missing")
image(fill1_svd, col=gray((0:100)/100), axes=FALSE, main="10 regressors")
image(fill2_svd, col=gray((0:100)/100), axes=FALSE, main="25 regressors")
image(fill3_svd, col=gray((0:100)/100), axes=FALSE, main="50 regressors")
par(opar_svd)
mtext("svd", side = 3, line = -21, outer = TRUE)

# It seems that the filling library will change the random seed, so it is set again.
set.seed(22200056)

# extract proteins data
q8.proteins <- global_data[, c(
    "kappa_casein",
    "alpha_s2_casein",
    "alpha_s1_casein",
    "beta_casein",
    "alpha_lactalbumin",
    "beta_lactoglobulin_a",
    "beta_lactoglobulin_b"
)]
# replace zeros with NA
for (col_idx in seq_len(ncol(q8.proteins))) {
    zero_inds <- which(q8.proteins[, col_idx] == 0.0)
    q8.proteins[zero_inds, col_idx] <- NA
}
# inpute matrix
q8.pca_inputed_proteins_1 <- impute(q8.proteins, method = "pca", k = 1, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_2 <- impute(q8.proteins, method = "pca", k = 2, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_3 <- impute(q8.proteins, method = "pca", k = 3, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_3 <- impute(q8.proteins, method = "pca", k = 4, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_3 <- impute(q8.proteins, method = "pca", k = 5, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_3 <- impute(q8.proteins, method = "pca", k = 6, tol = 1e-7, maxiter = 1000)
q8.pca_inputed_proteins_3 <- impute(q8.proteins, method = "pca", k = 7, tol = 1e-7, maxiter = 1000)

############################################################################
#                                Question 09                               #
############################################################################
q9.n_observations <- nrow(global_data)
# train set
q9.train_inds <- sample(
    1:q9.n_observations,
    size = 2 / 3 * q9.n_observations,
    replace = FALSE
)
q9.train <- data.frame(
    spectra = I(as.matrix(global_data[q9.train_inds, c(52:582)])),
    proteins = I(global_data[q9.train_inds, c(
        "kappa_casein",
        "alpha_s2_casein",
        "alpha_s1_casein",
        "beta_casein",
        "alpha_lactalbumin",
        "beta_lactoglobulin_a",
        "beta_lactoglobulin_b"
    )]),
    beta_ln_b = global_data[q9.train_inds, "beta_lactoglobulin_b"]
)
# test set, observations with zero beta_lactoglobulin_b values are removed
q9.test_inds <- c(1:q9.n_observations)[-q9.train_inds]
q9.test <- data.frame(
    spectra = I(as.matrix(global_data[q9.test_inds, c(52:582)])),
    beta_ln_b = global_data[q9.test_inds, "beta_lactoglobulin_b"]
)
q9.test <- q9.test[q9.test$beta_ln_b != 0, ]

# (a) stragety 1: omit observations with zero beta_lactoglobulin_b values
q9.train_omit <- q9.train[q9.train$beta_ln_b != 0.0, ]
# fitting
q9.pcr_omit <- pcr(
    formula = q9.train_omit$beta_ln_b ~ q9.train_omit$spectra,
    ncomp = 100,
    data = q9.train_omit,
    center = TRUE,
    scale = FALSE,
    validation = "CV"
)
print(summary(q9.pcr_omit))
plot(
    RMSEP(q9.pcr_omit),
    legendpos = "topleft",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
plot(q9.pcr_omit, ncomp = 45, asp = 1, line = TRUE)
# prediction
q9.rmse_omit <- rep(NA, 100)
for (i in 1:100) {
    q9.preds_omit <- predict(q9.pcr_omit, ncomp = i, newdata = q9.test$spectra)
    q9.rmse_omit[i] <- sqrt(mean((q9.preds_omit - q9.test$beta_ln_b)^2))
}
plot(
    q9.rmse_omit,
    type = "l",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
print(min(q9.rmse_omit))
print(q9.rmse_omit[45])

# (b) strategy 2: fill zero beta_lactoglobulin_b values with mean value
q9.train_mean <- q9.train
for (col_idx in seq_len(ncol(q9.train_mean$proteins))) {
    zero_inds <- which(q9.train_mean$proteins[, col_idx] == 0.0)
    q9.train_mean$proteins[zero_inds, col_idx] <- NA
}
q9.train_mean$proteins <- impute(q9.train_mean$proteins, method = "mean")
q9.train_mean$beta_ln_b <- q9.train_mean$proteins$beta_lactoglobulin_b
# fitting
q9.pcr_mean <- pcr(
    formula = q9.train_mean$beta_ln_b ~ q9.train_mean$spectra,
    ncomp = 100,
    data = q9.train_mean,
    center = TRUE,
    scale = FALSE,
    validation = "CV"
)
print(summary(q9.pcr_mean))
plot(
    RMSEP(q9.pcr_mean),
    legendpos = "topleft",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
plot(q9.pcr_mean, ncomp = 45, asp = 1, line = TRUE)
# prediction
q9.rmse_mean <- rep(NA, 100)
for (i in 1:100) {
    q9.preds_mean <- predict(q9.pcr_mean, ncomp = i, newdata = q9.test$spectra)
    q9.rmse_mean[i] <- sqrt(mean((q9.preds_mean - q9.test$beta_ln_b)^2))
}
plot(
    q9.rmse_mean,
    type = "l",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
print(min(q9.rmse_mean))
print(q9.rmse_mean[45])

# (c) strategy 3: fill zero beta_lactoglobulin_b values with pca imputed values
q9.train_pca <- q9.train
for (col_idx in seq_len(ncol(q9.train_pca$proteins))) {
    zero_inds <- which(q9.train_pca$proteins[, col_idx] == 0.0)
    q9.train_pca$proteins[zero_inds, col_idx] <- NA
}
q9.train_pca$proteins <- impute(q9.train_pca$proteins,
    method = "pca",
    k = 1,
    tol = 1e-7,
    maxiter = 1000
)
q9.train_pca$beta_ln_b <- q9.train_pca$proteins$beta_lactoglobulin_b
# fitting
q9.pcr_pca <- pcr(
    formula = q9.train_pca$beta_ln_b ~ q9.train_pca$spectra,
    ncomp = 100,
    data = q9.train_pca,
    center = TRUE,
    scale = FALSE,
    validation = "CV"
)
print(summary(q9.pcr_pca))
plot(
    RMSEP(q9.pcr_pca),
    legendpos = "topleft",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
plot(q9.pcr_pca, ncomp = 45, asp = 1, line = TRUE)
# prediction
q9.rmse_pca <- rep(NA, 100)
for (i in 1:100) {
    q9.preds_pca <- predict(q9.pcr_pca, ncomp = i, newdata = q9.test$spectra)
    q9.rmse_pca[i] <- sqrt(mean((q9.preds_pca - q9.test$beta_ln_b)^2))
}
plot(
    q9.rmse_pca,
    type = "l",
    xlab = "principal component",
    ylab = "Root Mean Squared Error (RMSE)"
)
print(min(q9.rmse_pca))
print(q9.rmse_pca[45])
