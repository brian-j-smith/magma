#################### Global Variables ####################

gpu <- TRUE

compare <- function(a, b, class="magma") {
   stopifnot(is(a, class), all.equal(dim(a), dim(b)),
             all.equal(as.numeric(a), as.numeric(b)))
}

#################### Tests: Constructors ####################

m <- 10
n <- 5
x <- rep(1, m)
A <- matrix(seq(m*n), m, n)

mx <- magma(x, gpu=gpu)
mA <- magma(A, gpu=gpu)

stopifnot(all.equal(magma(x, m, 1, gpu=gpu), mx))
stopifnot(all.equal(magma(mx, gpu=gpu), mx))

stopifnot(all.equal(as(mx, "numeric"), x))
stopifnot(all.equal(as(mx, "matrix"), as.matrix(x)))

stopifnot(all.equal(as(mA, "matrix"), A))
stopifnot(all.equal(as(mA, "numeric"), as.numeric(A)))

gpu(mx) <- !gpu
stopifnot(gpu(mx) == !gpu)
gpu(mx) <- gpu

compare(mA[1:2,], A[1:2,])
compare(mA[1:2,,drop=FALSE], A[1:2,,drop=FALSE])

compare(mA[,1:2], A[,1:2])
compare(mA[,1:2,drop=FALSE], A[,1:2,drop=FALSE])

compare(mA[,], A[,])
compare(mA[,,drop=FALSE], A[,,drop=FALSE])

compare(mA[1,], A[1,], class="numeric")
compare(mA[1,,drop=FALSE], A[1,,drop=FALSE])

compare(mA[,1], A[,1], class="numeric")
compare(mA[,1,drop=FALSE], A[,1,drop=FALSE])


#################### Tests: Addition ####################

m <- 2
n <- 3
x <- rep(1, m)
A <- matrix(seq(m*n), m, n)

mx <- magma(x, gpu=gpu)
mA <- magma(A, gpu=gpu)

y <- as.matrix(x + x)
compare(mx + mx, y)
compare(mx + x, y)
compare(x + mx, y)

B <- A + A
compare(mA + mA, B)
compare(mA + A, B)
compare(A + mA, B)


#################### Tests: Multiplication ####################

m <- 2
n <- 3
x <- rep(1, m)
y <- rep(1, n)
A <- matrix(seq(m*n), m, n)
B <- matrix(seq(m*n), n, m)

mx <- magma(x, gpu=gpu)
my <- magma(y, gpu=gpu)
mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)

C <- x %*% A;
## Conformable matrix operations
compare(t(mx) %*% mA, C)
compare(t(mx) %*% A, C)
compare(x %*% mA, C)

## Non-conformable operations
stopifnot(class(try(mx %*% mA, silent = TRUE)) == "try-error")
stopifnot(class(try(mx %*% A, silent = TRUE)) == "try-error")
stopifnot(class(try(mx %*% mB, silent = TRUE)) == "try-error")
stopifnot(class(try(mx %*% B, silent = TRUE)) == "try-error")
stopifnot(class(try(t(mx) %*% mB, silent = TRUE)) == "try-error")
stopifnot(class(try(t(mx) %*% B, silent = TRUE)) == "try-error")
stopifnot(class(try(t(my) %*% mA, silent = TRUE)) == "try-error")
stopifnot(class(try(t(my) %*% A, silent = TRUE)) == "try-error")
stopifnot(class(try(y %*% mA, silent = TRUE)) == "try-error")

C <- A %*% y;
## Conformable operations
compare(mA %*% my, C)
compare(mA %*% y, C)
compare(A %*% my, C)

## Non-conformable operations
stopifnot(class(try(mA %*% mx, silent = TRUE)) == "try-error")
stopifnot(class(try(mA %*% x, silent = TRUE)) == "try-error")
stopifnot(class(try(A %*% mx, silent = TRUE)) == "try-error")

C <- A %*% B;
## Conformable operations
compare(mA %*% mB, C)
compare(mA %*% B, C)
compare(A %*% mB, C)

## Non-conformable operations
stopifnot(class(try(mA %*% mA, silent = TRUE)) == "try-error")
stopifnot(class(try(mA %*% A, silent = TRUE)) == "try-error")
stopifnot(class(try(A %*% mA, silent = TRUE)) == "try-error")


#################### Tests: crossproduct ####################

m <- 2
n <- 3
x <- rep(1, m)
y <- rep(1, n)
A <- matrix(seq(m*n), m, n)
B <- matrix(seq(m*n), n, m)

mx <- magma(x, gpu=gpu)
my <- magma(y, gpu=gpu)
mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)

C <- crossprod(x, A)
## Conformable matrix operations
compare(crossprod(mx, mA), C)
compare(crossprod(mx, A), C)
compare(crossprod(x, mA), C)

## Non-conformable operations
stopifnot(class(try(crossprod(my, mA), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(my, A), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(y, mA), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(mx, mB), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(mx, B), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(x, mB), silent = TRUE)) == "try-error")

C <- crossprod(B, y)
## Conformable matrix operations
compare(crossprod(mB, my), C)
compare(crossprod(mB, y), C)
compare(crossprod(B, my), C)

## Non-conformable operations
stopifnot(class(try(crossprod(my, mA), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(my, A), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(y, mA), silent = TRUE)) == "try-error")

C <- crossprod(x, x)
## Conformable matrix operations
compare(crossprod(mx, mx), C)
compare(crossprod(mx, x), C)
compare(crossprod(x, mx), C)

## Non-conformable operations
stopifnot(class(try(crossprod(mx, my), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(mx, y), silent = TRUE)) == "try-error")
stopifnot(class(try(crossprod(x, my), silent = TRUE)) == "try-error")

C <- crossprod(A)
## Conformable matrix operations
compare(crossprod(mA), C)

C <- crossprod(x)
## Conformable matrix operations
compare(crossprod(mx), C)


#################### Tests: tcrossproduct ####################

m <- 2
n <- 3
x <- rep(1, m)
y <- rep(1, n)
A <- matrix(seq(m*n), m, n)
B <- matrix(seq(m*n), n, m)

mx <- magma(x, gpu=gpu)
my <- magma(y, gpu=gpu)
mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)

C <- tcrossprod(A, A)
## Conformable matrix operations
compare(tcrossprod(mA, mA), C)
compare(tcrossprod(mA, A), C)
compare(tcrossprod(A, mA), C)

## Non-conformable operations
stopifnot(class(try(tcrossprod(mB, mA), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(mB, A), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(B, mA), silent = TRUE)) == "try-error")

C <- tcrossprod(x, y)
## Conformable matrix operations
compare(tcrossprod(mx, my), C)
compare(tcrossprod(mx, y), C)
compare(tcrossprod(x, my), C)

## Non-conformable operations
stopifnot(class(try(tcrossprod(my, mA), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(my, A), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(mA, my), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(mA, y), silent = TRUE)) == "try-error")
stopifnot(class(try(tcrossprod(A, my), silent = TRUE)) == "try-error")

C <- tcrossprod(x, B)
## Conformable matrix operations
compare(tcrossprod(x, mB), C)

C <- tcrossprod(x, as.matrix(x))
## Conformable matrix operations
compare(tcrossprod(x, mx), C)

C <- tcrossprod(A)
## Conformable matrix operations
compare(tcrossprod(mA), C)

C <- tcrossprod(x)
## Conformable matrix operations
compare(tcrossprod(mx), C)


#################### Tests: General Linear Solver ####################

A <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
B <- matrix(c(1, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3)
y <- rep(1, nrow(A))

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)
my <- magma(y, gpu=gpu)

compare(solve(mA), solve(A))
compare(solve(mA, mB), solve(A, B))
compare(solve(mA, B), solve(A, B))
compare(solve(A, mB), solve(A, B))

compare(solve(mA, my), solve(A, as.matrix(y)))
compare(solve(mA, y), solve(A, y), class="numeric")
compare(solve(A, my), solve(A, my))

## Check large matrices
n <- 1000
x <- matrix(rnorm(n*n), n, n)
A <- tcrossprod(x)

mA <- magma(A, gpu=gpu)

compare(solve(mA), solve(A))

## Nonsingular and non-conformable matrices
m <- 2; n <- 3
A <- matrix(1, m, m)
B <- matrix(seq(m*n), n, m)

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)

stopifnot(class(try(solve(mA), silent = TRUE)) == "try-error")
stopifnot(class(try(solve(mB), silent = TRUE)) == "try-error")
stopifnot(class(try(solve(mA, mB), silent = TRUE)) == "try-error")
stopifnot(class(try(solve(mA, B), silent = TRUE)) == "try-error")
stopifnot(class(try(solve(A, mB), silent = TRUE)) == "try-error")


#################### Tests: Backward/Forward Solvers ####################

A <- matrix(c(1, 0, 0, 0.4, 1, 0, 0.2, 0.3, 1), 3, 3)
B <- matrix(c(1, 0, 0, 0.4, 1, 0, 0.2, 0.3, 1, 0.1, 0.1, 0.1), 3, 4)
I <- diag(nrow(A))
x <- rep(1, 3)

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)
mI <- magma(I, gpu=gpu)
mx <- magma(x, gpu=gpu)

compare(backsolve(mA, mI), backsolve(A, I))
compare(backsolve(mA, I), backsolve(A, I))
compare(backsolve(A, mI), backsolve(A, I))

compare(backsolve(mA, mx), backsolve(A, as.matrix(x)))
compare(backsolve(mA, x), backsolve(A, x), class="numeric")
compare(backsolve(A, mx), backsolve(A, as.matrix(x)))

compare(forwardsolve(t(mA), mI), forwardsolve(t(A), I))
compare(forwardsolve(t(mA), I), forwardsolve(t(A), I))
compare(forwardsolve(t(A), mI), forwardsolve(t(A), I))

compare(forwardsolve(t(mA), mx), forwardsolve(t(A), as.matrix(x)))
compare(forwardsolve(t(mA), x), forwardsolve(t(A), x), class="numeric")
compare(forwardsolve(t(A), mx), forwardsolve(t(A), as.matrix(x)))


#################### Tests: Cholesky Factorization ####################

A <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
B <- matrix(seq(9), 3, 3)

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)

compare(chol(mA), chol(A))
compare(chol2inv(chol(mA)), chol2inv(chol(A)))

stopifnot(class(try(chol(mB), silent = TRUE)) == "try-error")

## Check large matrices
n <- 1000
x <- matrix(rnorm(n*n), n, n)
A <- tcrossprod(x)

mA <- magma(A, gpu=gpu)

compare(chol(mA), chol(A))


#################### Tests: LU Factorization ####################

A <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
B <- matrix(c(1, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3)
y <- seq(nrow(A))

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)
my <- magma(y, gpu=gpu)

a <- lu(mA)

compare(solve(a), solve(A))
compare(solve(a, y), solve(A, y), class="numeric")
compare(solve(a, B), solve(A, B))


#################### Tests: QR Factorization ####################

A <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), 3, 3)
B <- matrix(c(1, 0, 0, 0, 0, 1, 0, 1, 0), 3, 3)
y <- seq(nrow(A))

mA <- magma(A, gpu=gpu)
mB <- magma(B, gpu=gpu)
my <- magma(y, gpu=gpu)

a <- qr(mA)
b <- qr(A, LAPACK=TRUE)

compare(solve(a), solve(b))
compare(solve(a, y), solve(b, y), class="numeric")
compare(solve(a, B), solve(b, B))

