# function to use in both estimate_mupp_params and estimate_mupp_thetas
#' @importFrom stats
#'             setNames
#'             as.formula
#' @importFrom utils
#'             head tail
estimate_mupp_header_ <- function(resp,
                                  items,
                                  type = "person"){

  # argument checks #

  # converting to lowercase
  resp  %<>% setNames(tolower(names(.)))
  items %<>% setNames(tolower(names(.)))

  # determine template function based in estimation algorithm
  if(type == "person"){
    template_fun <- function(...){
      c(simulate_mupp_resp(...)["resp"], list(...))
    } # END template_fun FUNCTION
  } else{
    template_fun <- simulate_mupp_resp
  } # END ifelse STATEMENT

  # pulling out template resp/item
  template <- lapply(do.call(what = template_fun,
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
  cast_resp      <- as.formula(paste0(person_name, "~", item_name))

  resp_adj       <- dcast(data      = as.data.table(resp_adj),
                          formula   = cast_resp,
                          value.var = resp_name)[ , -1] %>%
                    as.matrix()

  # updating to items/params/resp (if required)
  if(type == "person"){
    params_adj    <- as.matrix(items_adj[4:ncol(items_adj)])
    items_adj     <- as.matrix(items_adj[seq_len(3)])
  } else{
    params_adj    <- NULL
  } # END ifelse STATEMENT

  # return
  return(list(resp       = resp,
              items      = items,
              resp_adj   = resp_adj,
              items_adj  = items_adj,
              params_adj = params_adj))
} # END estimate_mupp_header_ FUNCTION
