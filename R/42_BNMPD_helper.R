load_model = function(from_env) {
  # if (is.null(to_env)) to_env <- new.env()
  to_env <- parent.frame()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
