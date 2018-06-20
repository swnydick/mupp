#' Estimate MUPP Parameters
#'
#' Estimate MUPP statement and person parameters given item responses and item
#' properties using MCMC
#'
#' @param resp a data.frame of (at least) [person, item, resp]
#' @param items a data.frame of (at least) [item, statement, dim]
#' @param method the estimation method (MCMC is the only one that works now)
#' @param ... other parameters to add later
#'
#' @return a list of [theta, params] arrays with the third dimension indicating
#'         the iteration
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames
#' @importFrom data.table dcast as.data.table
#' @export
estimate_mupp_params <- function(resp,
                                 items,
                                 method = "MCMC",
                                 ...){

  # argument checks #

  # converting to lowercase
  resp  %<>% setNames(tolower(names(.)))
  items %<>% setNames(tolower(names(.)))

  # pulling out template resp/item
  template <- lapply(do.call(what = simulate_mupp_resp,
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

  # run algorithm
  if(any(tolower(method) == "mcmc")){
    out <- estimate_mupp_params_mcmc(resp  = resp_adj,
                                     items = items_adj,
                                     ...)
  } else{
    stop(method, " method not implemented at this time.")
  } # END ifelse STATEMENT

  return(out)

} # END estimate_mupp_params FUNCTION
