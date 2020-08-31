# Tests ResidualMatrix multiplication.
# library(testthat); library(ResidualMatrix); source("setup.R"); source("test-mult.R")

expect_equal_product <- function(x, y) {
    expect_s4_class(x, "DelayedMatrix")
    X <- as.matrix(x)

    # standardize NULL dimnames.
    if (all(lengths(dimnames(X))==0L)) dimnames(X) <- NULL
    if (all(lengths(dimnames(y))==0L)) dimnames(y) <- NULL
    expect_equal(X, y)
}

test_that("ResidualMatrix right multiplication works as expected", {
    possibles <- spawn_scenarios(100, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), ncol=10)
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)
    }
})

test_that("ResidualMatrix left multiplication works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), nrow=10)
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)
    }
})

test_that("ResidualMatrix dual multiplication works as expected", {
    possibles1 <- spawn_scenarios(10, 20)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(test1$res %*% test2$res, test1$ref %*% test2$ref)

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(test1$res[0,] %*% test2$res, test1$ref[0,] %*% test2$ref)
            expect_equal_product(test1$res %*% test2$res[,0], test1$ref %*% test2$ref[,0])
            expect_equal_product(test1$res[,0] %*% test2$res[0,], test1$ref[,0] %*% test2$ref[0,])
            expect_equal_product(test1$res[0,] %*% test2$res[,0], test1$ref[0,] %*% test2$ref[,0])
        }
    }
})

##########################

test_that("ResidualMatrix lonely crossproduct works as expected", {
    possibles <- spawn_scenarios(90, 30)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res
        expect_equal_product(crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("ResidualMatrix crossproduct from right works as expected", {
    possibles <- spawn_scenarios(60, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))
    }
})

test_that("ResidualMatrix crossproduct from left works as expected", {
    possibles <- spawn_scenarios(40, 100)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))
    }
})

test_that("ResidualMatrix dual crossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(crossprod(test1$res, test2$res), crossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(crossprod(test1$res[,0], test2$res), crossprod(test1$ref[,0], test2$ref))
            expect_equal_product(crossprod(test1$res, test2$res[,0]), crossprod(test1$ref, test2$ref[,0]))
            expect_equal_product(crossprod(test1$res[0,], test2$res[0,]), crossprod(test1$ref[0,], test2$ref[0,]))
        }
    }
})

##########################

test_that("ResidualMatrix lonely tcrossproduct works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res
        expect_equal_product(tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("ResidualMatrix tcrossproduct from right works as expected", {
    possibles <- spawn_scenarios(60, 70)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector (this doesn't work).
        z <- rnorm(ncol(ref.y))
        expect_error(tcrossprod(bs.y, z), "non-conformable")
        expect_error(tcrossprod(ref.y, z), "non-conformable")

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))
    }
})

test_that("ResidualMatrix tcrossproduct from left works as expected", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))
    }
})

test_that("ResidualMatrix dual tcrossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(25, 50)
        for (test2 in possibles2) {

            expect_equal_product(tcrossprod(test1$res, test2$res), tcrossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(tcrossprod(test1$res[0,], test2$res), tcrossprod(test1$ref[0,], test2$ref))
            expect_equal_product(tcrossprod(test1$res, test2$res[0,]), tcrossprod(test1$ref, test2$ref[0,]))
            expect_equal_product(tcrossprod(test1$res[,0], test2$res[,0]), tcrossprod(test1$ref[,0], test2$ref[,0]))
        }
    }
})
