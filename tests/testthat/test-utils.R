test_that("%||% returns x when not NULL, y otherwise", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
})

test_that("ev_abort emits cli_abort with class", {
  expect_error(
    ev_abort("test message", class = "ev_test_class"),
    class = "ev_test_class"
  )
})

test_that("ev_inform emits cli_inform with class", {
  expect_message(
    ev_inform("test message", class = "ev_test_class"),
    regexp = "test message"
  )
})
