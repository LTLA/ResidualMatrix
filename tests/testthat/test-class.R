# Tests ResidualMatrix logic.
# library(testthat); library(ResidualMatrix); source("setup.R"); source("test-class.R")

set.seed(100002)
test_that("ResidualMatrix theory behaves as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_s4_class(test$res, "ResidualMatrix")
        expect_identical(test$res, ResidualMatrix(DelayedArray::seed(test$res)))

        expect_identical(dim(test$res), dim(test$ref))
        expect_equal(extract_array(test$res, list(1:10, 1:10)), test$ref[1:10, 1:10])
        expect_equal(extract_array(test$res, list(1:10, NULL)), test$ref[1:10,])
        expect_equal(extract_array(test$res, list(NULL, 1:10)), test$ref[,1:10])
        expect_equal(purgenames(as.matrix(test$res)), purgenames(test$ref))
    }
})

set.seed(100002)
test_that("ResidualMatrix theory behaves as expected with kept coefficients", {
    NR <- 1501
    NC <- 121

    # Run through a host of different design matrices.
    for (it in 2:5) {
        y <- matrix(rnorm(NR*NC), ncol=NC)
        design <- generate_design(nrow(y), it)

        for (k in 1:ncol(design)) {
            res <- ResidualMatrix(y, design, keep=k)

            fit <- lm.fit(x=design, y=as.matrix(y))
            ref <- design[,k,drop=FALSE] %*% fit$coefficients[k,,drop=FALSE] + fit$residuals
            dimnames(ref) <- NULL

            expect_equal(as.matrix(res), ref)
        }

        res <- ResidualMatrix(y, design, keep=seq_len(ncol(design)))
        expect_equal(as.matrix(res), y)
    }
})

set.seed(100002)
test_that("ResidualMatrix theory behaves as expected with restriction", {
    NR <- 501
    NC <- 121
    restrict <- sample(NR, NR/2)

    # Run through a host of different design matrices.
    for (it in 2:5) {
        y <- matrix(rnorm(NR*NC), ncol=NC)
        design <- generate_design(nrow(y), it)

        res <- ResidualMatrix(y, design, restrict=restrict)
        ref <- ResidualMatrix(y[restrict,,drop=FALSE], design[restrict,,drop=FALSE])
        expect_equal(as.matrix(res[restrict,,drop=FALSE]), as.matrix(ref))

        fit <- lm.fit(x=design[restrict,,drop=FALSE], y=as.matrix(y[restrict,,drop=FALSE]))
        ref <- as.matrix(y) - design %*% fit$coefficients
        dimnames(ref) <- NULL

        expect_equal(as.matrix(res), ref)
        res <- ResidualMatrix(y, design, restrict=!logical(NR))
        expect_equal(as.matrix(res), as.matrix(ResidualMatrix(y, design)))

        # Works correctly when combined with keep.
        for (k in 1:ncol(design)) {
            res <- ResidualMatrix(y, design, keep=k, restrict=restrict)

            test <- fit$coefficients
            test[k,] <- 0
            ref <- as.matrix(y) - design %*% test
            dimnames(ref) <- NULL

            expect_equal(as.matrix(res), ref)
        }
    }
})

set.seed(100002)
test_that("ResidualMatrix theory throws with non-full rank", {
    y <- matrix(rnorm(1000), nrow=20)
    g <- gl(4, 5)
    g2 <- rep(gl(2, 2), 5)
    X <- cbind(model.matrix(~g2), model.matrix(~g))
    expect_error(res <- ResidualMatrix(y, X), "full rank")
})
