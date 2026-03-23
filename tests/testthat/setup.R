op <- options(quickcheck.tests = 20L, quickcheck.shrinks = 20L)

withr::defer(options(op), teardown_env())
