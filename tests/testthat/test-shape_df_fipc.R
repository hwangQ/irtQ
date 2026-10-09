# shape_df_fipc() builds the metadata of fixed and new items for FIPC

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_fix <- bring.flexmirt(file = prm_file, "par")$Group1$full_df[1:10, ]

test_that("shape_df_fipc() repeats a single cats and model for every new item", {
  new_ids <- paste0("N", 1:4)
  # single values
  meta1 <- shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = 2, model = "3PLM")
  # the same input written out in full
  meta2 <- shape_df_fipc(
    x = x_fix, fix.loc = 1:10, item.id = new_ids,
    cats = rep(2, 4), model = rep("3PLM", 4)
  )
  expect_identical(meta1, meta2)
  expect_equal(nrow(meta1), 14L)
  expect_identical(as.character(meta1$id[11:14]), new_ids)
  expect_true(all(meta1$model[11:14] == "3PLM"))
})

test_that("shape_df_fipc() stops when cats or model do not match the number of new items", {
  new_ids <- paste0("N", 1:3)
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = c(2, 2), model = "3PLM"),
    "must be 1 or equal to the number of new items"
  )
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = 1:10, item.id = new_ids, cats = 2, model = c("3PLM", "2PLM")),
    "must be 1 or equal to the number of new items"
  )
})

test_that("shape_df_fipc() stops when fix.loc is missing or invalid", {
  new_ids <- paste0("N", 1:4)
  # fix.loc is not given
  expect_error(
    shape_df_fipc(x = x_fix, item.id = new_ids, cats = 2, model = "3PLM"),
    "fix.loc"
  )
  # fix.loc is shorter than the number of fixed items
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = 1:3, item.id = new_ids, cats = 2, model = "3PLM"),
    "fix.loc"
  )
  # fix.loc repeats a position
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = c(1:9, 9), item.id = new_ids, cats = 2, model = "3PLM"),
    "fix.loc"
  )
  # fix.loc points beyond the final form
  expect_error(
    shape_df_fipc(x = x_fix, fix.loc = c(1:9, 15), item.id = new_ids, cats = 2, model = "3PLM"),
    "fix.loc"
  )
})

test_that("shape_df_fipc() places the fixed and new items at their positions", {
  meta <- shape_df_fipc(
    x = x_fix[1:3, ], fix.loc = c(1, 3, 5), item.id = c("N1", "N2"),
    cats = 2, model = "3PLM"
  )
  expect_identical(as.character(meta$id), c(x_fix$id[1], "N1", x_fix$id[2], "N2", x_fix$id[3]))
})
