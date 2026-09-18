.rgse_net_connect_timeout <- 60
.rgse_net_stall_bytes     <- 1
.rgse_net_stall_seconds   <- 180
.rgse_net_total_timeout   <- 1800
.rgse_net_base_timeout    <- 600
.rgse_net_attempts        <- 4

.rgse_net_configure <- function(verbose = TRUE) {
  options(timeout = .rgse_net_base_timeout)

  options(RCurlOptions = utils::modifyList(
    if (is.list(getOption("RCurlOptions"))) getOption("RCurlOptions") else list(),
    list(connecttimeout = .rgse_net_connect_timeout,
         timeout        = .rgse_net_total_timeout,
         low.speed.limit = .rgse_net_stall_bytes,
         low.speed.time  = .rgse_net_stall_seconds)))

  if (requireNamespace("httr", quietly = TRUE)) {
    tryCatch(
      httr::set_config(httr::config(connecttimeout  = .rgse_net_connect_timeout,
                                    timeout         = .rgse_net_total_timeout,
                                    low_speed_limit = .rgse_net_stall_bytes,
                                    low_speed_time  = .rgse_net_stall_seconds),
                       override = FALSE),
      error = function(e) invisible(NULL))
  }

  if (verbose)
    print(paste0("Network guards: connect ", .rgse_net_connect_timeout, "s, abort after ",
                 .rgse_net_stall_seconds, "s below ", .rgse_net_stall_bytes,
                 " byte/s, hard cap ", .rgse_net_total_timeout, "s, ",
                 .rgse_net_attempts, " attempts"))
  invisible(NULL)
}

.rgse_net_retry <- function(expr, attempts = .rgse_net_attempts, label = NULL, quiet = FALSE) {
  call_expr <- substitute(expr)
  env <- parent.frame()
  for (k in seq_len(attempts)) {
    res <- tryCatch(eval(call_expr, env), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (k == attempts) {
      if (!quiet)
        print(paste0("Network call failed after ", attempts, " attempts",
                     if (!is.null(label)) paste0(" (", label, ")") else "", ": ",
                     conditionMessage(res)))
      return(NULL)
    }
    wait <- min(120, 5 * 2^(k - 1)) + stats::runif(1, 0, 5)
    if (!quiet)
      print(paste0("Network call failed",
                   if (!is.null(label)) paste0(" (", label, ")") else "",
                   ", retrying in ", round(wait), "s [", k, "/", attempts - 1, "]: ",
                   conditionMessage(res)))
    Sys.sleep(wait)
  }
  NULL
}
