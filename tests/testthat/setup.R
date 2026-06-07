# Limit BLAS/OpenMP thread usage during testing
#
# irtQ relies heavily on Rfast and Matrix for linear algebra operations
# (e.g., EM.R, est_irt.R, info.R, covirt.R). On some CRAN check machines,
# the underlying BLAS/LAPACK library (e.g., OpenBLAS) automatically spawns
# multiple threads for matrix computations. This causes the reported CPU
# time to be several times larger than the elapsed (wall-clock) time during
# `R CMD check` (e.g., "Running R code in 'testthat.R' had CPU time 6.7
# times elapsed time"), which CRAN flags as a NOTE because package checks
# should not make heavy use of multiple cores/threads.
#
# Forcing all relevant thread-control environment variables to "1" here
# ensures BLAS/OpenMP operations run single-threaded for the duration of
# the test suite, without changing any package code or user-facing behavior.
Sys.setenv(OMP_NUM_THREADS = "1")        # OpenMP thread limit
Sys.setenv(OMP_THREAD_LIMIT = "1")       # OpenMP thread limit (alternate variable)
Sys.setenv(OPENBLAS_NUM_THREADS = "1")   # OpenBLAS thread limit
Sys.setenv(MKL_NUM_THREADS = "1")        # Intel MKL thread limit
Sys.setenv(BLIS_NUM_THREADS = "1")       # BLIS thread limit
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1") # Apple Accelerate/vecLib thread limit
