#' Estimate MUPP Thetas
#'
#' Estimate MUPP person parameters given item responses, item parameters, and
#' item properties.
#'
#' @param resp a data.frame of (at least) [person, item, resp]
#' @param items a data.frame of (at least) [item, statement, dim, alpha, delta, tau]
#' @param method the estimation method ("newton_r", "bfgs_r", "newton_c", "bfgs_c",
#'        "MCMC")
#' @param ... other parameters to pass to the method
#'
#' @return a list of [thetas, sems, hessians]
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames
#' @importFrom data.table dcast as.data.table
#' @export
estimate_mupp_thetas <- function(resp,
                                 items,
                                 method  = c("newton_r", "bfgs_r",
                                             "newton_c", "bfgs_c",
                                             "MCMC"),
                                 control = list(),
                                 ...){

  # argument checks #

  # fix method
  method   <- match.arg(method)

  # converting to lowercase
  resp   %<>% setNames(tolower(names(.)))
  items  %<>% setNames(tolower(names(.)))

  # pulling out template resp/item
  template <- lapply(do.call(what = function(...){
                               c(simulate_mupp_resp(...)["resp"], list(...))
                             },
                             args = simulate_mupp_params()),
                     FUN = names)

  # ordering the columns
  resp   %<>% check_names(template$resp)
  items  %<>% check_names(template$items)

  # extracting variable names

  # overall names
  resp_names     <- names(resp)
  item_names     <- names(items)

  # individual names
  person_name    <- head(resp_names, 1)
  item_name      <- head(item_names, 1)
  statement_name <- item_names[2]
  resp_name      <- tail(resp_names, 1)

  # pulling out old item/statement/person and setting new item to index
  items_adj      <- sequence_column(df     = items,
                                    column = item_name) %>%
                    sequence_column(column = statement_name)
  resp_adj       <- sequence_column(df     = resp,
                                    column = person_name) %>%
                    sequence_column(column = item_name,
                                    old_values = unique(items[[item_name]]))

  # transforming response to be of the appropriate form
  cast_resp     <- as.formula(paste0(person_name, "~", item_name))

  resp_adj      <- dcast(data      = as.data.table(resp_adj),
                         formula   = cast_resp,
                         value.var = resp_name)[ , -1] %>%
                   as.matrix()

  # updating to items/params/resp
  params_adj    <- as.matrix(items_adj[4:ncol(items_adj)])
  items_adj     <- as.matrix(items_adj[seq_len(3)])

  # run algorithm
  algorithm_fun <- switch(method,
                          newton_r = estimate_mupp_thetas_newton_r,
                          stop(method, " method not implemented at this time."))

  out           <- algorithm_fun(resp    = resp_adj,
                                 params  = params_adj,
                                 items   = items_adj,
                                 control = control,
                                 ...)

  # calculate SEM and return

  # Note: most of this is similar to estimate_mupp_params ... maybe combine similar parts?

  return(out)

} # END estimate_mupp_thetas FUNCTION
