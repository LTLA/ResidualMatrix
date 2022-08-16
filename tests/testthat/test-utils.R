# Tests ResidualMatrix utilities.
# library(testthat); library(ResidualMatrix); source("setup.R"); source("test-utils.R")

set.seed(100001)
test_that("ResidualMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_equal(rowSums(test$res), rowSums(test$ref))
        expect_equal(colSums(test$res), colSums(test$ref))
        expect_equal(rowMeans(test$res), rowMeans(test$ref))
        expect_equal(colMeans(test$res), colMeans(test$ref))

        tdef <- t(test$res)
        expect_s4_class(tdef, "ResidualMatrix") # still a ResMat!
        expect_equal(t(tdef), test$res)
        expect_equal(purgenames(as.matrix(tdef)), purgenames(t(test$ref)))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$res)))
        colnames(test$res) <- spawn_names
        expect_identical(spawn_names, colnames(test$res))
        expect_s4_class(test$res, "ResidualMatrix") # still a ResMat!
    }
})

set.seed(10000101)
test_that("ResidualMatrix silly inputs work as expected", {
    default <- ResidualMatrix()
    expect_identical(dim(default), c(0L, 0L))
    val <- as.matrix(default)
    dimnames(val) <- NULL
    expect_identical(val, matrix(0,0,0))

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(ResidualMatrix(y, design=cbind(1)), "nrow.* should be equal")
})

set.seed(1000011)
test_that("ResidualMatrix subsetting works as expected", {
    expect_equal_and_resmat <- function(x, y) {
        expect_s4_class(x, "ResidualMatrix") # class is correctly preserved by direct seed modification.
        expect_equal(purgenames(as.matrix(x)), purgenames(y))
    }

    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$res))
        j <- sample(ncol(test$res))
        expect_equal_and_resmat(test$res[i,], test$ref[i,])
        expect_equal_and_resmat(test$res[,j], test$ref[,j])
        expect_equal_and_resmat(test$res[i,j], test$ref[i,j])

        # Works with zero dimensions.
        expect_equal_and_resmat(test$res[0,], test$ref[0,])
        expect_equal_and_resmat(test$res[,0], test$ref[,0])
        expect_equal_and_resmat(test$res[0,0], test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_equal(test$res[i[1],], test$ref[i[1],])
        expect_equal(test$res[,j[1]], test$ref[,j[1]])
        expect_equal_and_resmat(test$res[i[1],drop=FALSE], test$ref[i[1],,drop=FALSE])
        expect_equal_and_resmat(test$res[,j[1],drop=FALSE], test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$res)
        expect_equal(t(alt[,i]), test$res[i,])
        expect_equal(t(alt[j,]), test$res[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$res)))
        colnames(test$res) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_equal_and_resmat(test$res[,ch], test$ref[,ch])
    }
})
