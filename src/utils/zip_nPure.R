zip_nPure = function(.x, .fields = NULL, .simplify = FALSE) {
  if (length(.x) == 0)
    return(list())
  if (is.null(.fields)) {
    if (is.null(names(.x[[1]]))) {
      .fields <- seq_along(.x[[1]])
    } else {
      .fields <- stats::setNames(names(.x[[1]]), names(.x[[1]]))
    }
  } else {
    if (is.character(.fields) && is.null(names(.fields))) {
      names(.fields) <- .fields
    }
  }
  out <- lapply(.fields, function(i) lapply(.x, .subset2, i))
  if (.simplify)
    out <- lapply(out, simplify_if_possible)
  out
}