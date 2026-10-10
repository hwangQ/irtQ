# panel_info() extracts the stages and pathways of an MST route map

test_that("panel_info() returns the stages and pathways of the simMST 1-3-3 panel", {
  pi <- panel_info(simMST$route_map)
  expect_equal(pi$config, list(stage.1 = 1L, stage.2 = 2:4, stage.3 = 5:7))
  expect_equal(unname(pi$pathway),
               matrix(c(1, 2, 5,  1, 2, 6,  1, 3, 5,  1, 3, 6,  1, 3, 7,  1, 4, 6,  1, 4, 7),
                      ncol = 3, byrow = TRUE))
  expect_equal(unname(pi$n.module), c(1L, 3L, 3L))
  expect_equal(pi$n.stage, 3L)
  # a data frame gives the same result
  expect_equal(panel_info(as.data.frame(simMST$route_map)), pi)
})

test_that("panel_info() finds every stage when module indices do not follow the stages", {
  # stage 2 modules are numbered 5 to 7 and stage 3 modules 2 to 4
  rm <- matrix(0, 7, 7)
  rm[1, 5:7] <- 1
  rm[5, 2:3] <- 1
  rm[6, 2:4] <- 1
  rm[7, 3:4] <- 1
  pi <- panel_info(rm)
  expect_equal(pi$config, list(stage.1 = 1L, stage.2 = 5:7, stage.3 = 2:4))
  expect_equal(nrow(pi$pathway), 7L)
})

test_that("panel_info() stops on route maps that do not define an MST panel", {
  # a module without any transition (module 4) leaves a stage without a next stage
  rm <- matrix(0, 4, 4)
  rm[1, 2:3] <- 1
  expect_error(panel_info(rm), "next stage")
  # a cycle
  cyc <- matrix(0, 3, 3)
  cyc[1, 2] <- 1
  cyc[2, 3] <- 1
  cyc[3, 2] <- 1
  expect_error(panel_info(cyc), "cycle")
  # not binary or not square
  expect_error(panel_info(matrix(2, 2, 2)), "binary")
  expect_error(panel_info(matrix(0, 2, 3)), "square")
  expect_error(panel_info(matrix(0, 1, 1)), "two stages")
})
