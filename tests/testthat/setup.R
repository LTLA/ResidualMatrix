generate_design <- function(nobs, code) {
    if (code==1L) {
        design <- NULL
    } else if (code==2L) {
        cov <- rnorm(nobs)
        design <- model.matrix(~cov)
    } else if (code==3L) {
        g <- factor(rep(1:3, length.out=nobs))
        design <- model.matrix(~0 + g)
    } else if (code==4L) {
        cov <- rnorm(nobs)
        g <- factor(rep(1:2, length.out=nobs))
        design <- model.matrix(~0 + g + cov)
    } else if (code==5L) {
        design <- cbind(rnorm(nobs))
    }
    design
}

spawn_scenarios_basic <- function(NR, NC, CREATOR, REALIZER) {
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:5) {
            if (trans) {
                # Ensure output matrix has NR rows and NC columns after t().
                y <- CREATOR(NC, NR)
            } else {
                y <- CREATOR(NR, NC)
            }
            ref <- REALIZER(y)

            # Run through a host of different design matrices.
            design <- generate_design(nrow(y), it)

            res <- ResidualMatrix(y, design)
            if (is.null(design)) {
                ref <- as.matrix(y)
                ref <- sweep(ref, 2, colMeans(ref), "-")
            } else {
                ref <- lm.fit(x=design, y=as.matrix(y))$residuals
            }

            if (trans) {
                res <- t(res)
                ref <- t(ref)
            }

            output[[counter]] <- list(res=res, ref=ref)
            counter <- counter+1L
        }
    }
    output
}

spawn_scenarios <- function(NR=50, NC=20) {
    c(
        spawn_scenarios_basic(NR, NC,
            CREATOR=function(r, c) {
                matrix(rnorm(r*c), ncol=c)
            },
            REALIZER=identity
        ),
        spawn_scenarios_basic(NR, NC,
            CREATOR=function(r, c) {
                Matrix::rsparsematrix(r, c, 0.1)
            },
            REALIZER=as.matrix
        )
    )
}

purgenames <- function(mat) {
    if (identical(dimnames(mat), list(NULL, NULL))) {
        dimnames(mat) <- NULL
    }
    mat
}
