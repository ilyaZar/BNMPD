copy_env = function(to_env = NULL, from_env) {
  if (is.null(to_env)) to_env <- new.env()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
copy_env_to_parent = function(from_env) {
  to_env <- parent.frame()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
