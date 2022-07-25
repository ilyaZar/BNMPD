#' #' @description Helper function to retrieve US state
#' #'   names to be used.
#' #'
#' #' @details First, states that have a-typical
#' #'   behavior are excluded: US (which is the overall
#' #'   USA so actually no states), then states that
#' #'   have permanent zeros (never produce a particular
#' #'   type of energy), then 'NJ' (which is just weird)
#' #'   and finally CO (just like that, no particular
#' #'   reason). The remaining states are retrieved in
#' #'   order of their population, see
#' #'   \code{num_biggest} for details.
#' #'
#' #' @param num_biggest numeric index vector specifying
#' #'   the indices (first to last ) number of states
#' #'   names to retrieve, starting with the most
#' #'   populated. So \code{num_biggest = 1:10} returns
#' #'   the names of the first 10 biggest US states
#' #'   (according to population).
#' #'
#' #' @return a character vector of states that
#' #'   satisfies above conditions.
#' get_state_names = function(num_biggest) {
#'
#'   # states_all1 <- unique(data$state)
#'   # states_all1 <- states_all1[states_all1 != "US"]
#'   states_all2 <- c("Alabama",
#'                    #"Alaska",
#'                    "Arizona","Arkansas","California",
#'                    "Colorado","Connecticut","Delaware","District of Columbia",
#'                    "Florida","Georgia", #"Hawaii",
#'                    "Idaho",
#'                    "Illinois","Indiana","Iowa","Kansas", "Kentucky","Louisiana",
#'                    "Maine","Maryland","Massachusetts", "Michigan","Minnesota",
#'                    "Mississippi","Missouri","Montana", "Nebraska","Nevada",
#'                    "New Hampshire","New Jersey","New Mexico", "New York",
#'                    "North Carolina","North Dakota","Ohio","Oklahoma", "Oregon",
#'                    "Pennsylvania","Rhode Island","South Carolina",
#'                    "South Dakota","Tennessee","Texas","Utah","Vermont","Virginia",
#'                    "Washington","West Virginia","Wisconsin","Wyoming")
#'   state_codes <- c("AL",
#'                    #"AK",
#'                    "AZ","AR","CA","CO","CT", "DE","DC","FL","GA",
#'                    #"HI",
#'                    "ID", "IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN",
#'                    "MS", "MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH",
#'                    "OK", "OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA",
#'                    "WV", "WI", "WY")
#'   names(state_codes) <- states_all2
#'   states_ordered <- c("California","Texas","Florida","New York","Pennsylvania",
#'                       "Illinois","Ohio",  "Georgia","North Carolina","Michigan",
#'                       "New Jersey","Virginia","Washington", "Arizona",
#'                       "Massachusetts","Tennessee","Indiana","Missouri",
#'                       "Maryland", "Wisconsin","Colorado","Minnesota",
#'                       "South Carolina","Alabama","Louisiana", "Kentucky",
#'                       "Oregon","Oklahoma","Connecticut","Utah","Iowa","Nevada",
#'                       "Arkansas","Mississippi","Kansas","New Mexico","Nebraska",
#'                       "West Virginia", "Idaho",
#'                       #"Hawaii",
#'                       "New Hampshire","Maine",
#'                       "Montana","Rhode Island","Delaware","South Dakota",
#'                       "North Dakota",
#'                       #"Alaska",
#'                       "District of Columbia","Vermont",
#'                       "Wyoming")
#'   states_to_use  <- state_codes[states_ordered]
#'
#'   states_exlucde <- c("AK", "DC", "DE", "HI", "ID", "IN", "KY",
#'                       "MS", "MT", "ND", "NM", "NV", "OK", "RI",
#'                       "SD", "UT", "WV", "WY",
#'                       "US",
#'                       "NJ",
#'                       "CO")
#'   states_to_use <- setdiff(states_to_use, states_exlucde)[1:num_biggest]
#'   return(states_to_use)
#' }

#'
#'
#' #' @description Returns initialization for latent
#' #'   states
#' #'
#' #' @details This is required for the SMC/particle
#' #'   filter; the current method infers starting
#' #'   values for the states from the data: zeros are
#' #'   set to a small value termed 'zero_lower_bound',
#' #'   while the other state components are
#' #'   proportional to the log of the scaled
#' #'   measurement y_t of the corresponding energy
#' #'   fraction (the \code{d} component). Scaling of
#' #'   the \code{y_t[t, d, n]} is done with a scaling
#' #'   factor before taking logs (by default set to 1).
#' #'
#' #' @param y_t the measurements e.g. for the dirichlet multinomial the raw energy
#' #'   counts
#' #' @param NN cross sectional dimension
#' #' @param TT time series dimension
#' #' @param DD number of components in the measurement vector
#' #'   (dimension of \code{y_t} i.e. number of energy types)
#' #'
#' #' @return an array of dimension \code{TTxDDxNN} containing the initialization
#' #'   values for the latent states
#' get_states_init = function(y_t, NN, TT, DD) {
#'   zero_lower_bound <- 0.001
#'   state_scale <- rep(1, times = DD)
#'   states_init <- array(0, c(TT, DD, NN))
#'   options(warn = 2)
#'   for (i in 1:NN) {
#'     # if (i == 11) browser()
#'     for (d in 1:DD) {
#'       states_init_temp <- abs(y_t[, d, i])
#'       states_init_temp[states_init_temp == 0] <- zero_lower_bound
#'       # print(i)
#'       # print(d)
#'       states_init_temp <- tryCatch(log(states_init_temp/state_scale[d]))
#'       states_init[, d, i] <- states_init_temp
#'     }
#'   }
#'   options(warn = 0)
#'   return(states_init)
#' }
#' )
#' #' @description Retrieves necessary subset of the
#' #'   data to be used for subsequent estimations.
#' #'
#' #' @param data_set full data set (as a data frame or
#' #'   tibble)
#' #' @param Y_to_use dependent variable (names and
#' #'   their labels) to use
#' #' @param Z_to_use a list of dimension equal to the
#' #'   number of dirichlet shares/fractions i.e. with
#' #'   each component being a vector of regressor names
#' #'   for the Z-type regressors (regressors for which
#' #'   parameters/coefficients are shared among all
#' #'   US-states)
#' #' @param U_to_use a list of dimension equal to the
#' #'   number of dirichlet shares/fractions i.e. with
#' #'   each component being a vector of regressor names
#' #'   for the U-type regressors (regressors for which
#' #'   parameters/coefficients vary among all
#' #'   US-states)
#' #'
#' #' @return a list with 6 components:
#' #'   \itemize{
#' #'   \item{y_t}{: an array of dimension
#' #'   \code{TTxDDxNN} (time series length, number of
#' #'   components, and number of cross sectional units)
#' #'   containing the measurements e.g. the raw energy
#' #'   counts or dirichlet shares}
#' #'   \item{num_counts:}{a
#' #'   matrix of dimension \code{TTxNN}; the second
#' #'   part of the measurements (for some distributions
#' #'   necessary e.g. the dirichlet multinomial) e.g.
#' #'   the number of total energy counts per time and
#' #'   per cross sectional unit}
#' #'   \item{Z:}{ an array of dimension
#' #'   \code{TTx(DD*num_regressors_per_component)xNN}
#' #'   containing the Z-type regressors; for each
#' #'   component \code{DD}; there might be a different
#' #'   number of regressors used so the number of
#' #'   columns has the component
#' #'   \code{num_regressors_per_component}; there are
#' #'   \code{NN} different cross sectional units (US
#' #'   states), so the final array dimension is
#' #'   \code{NN}}
#' #'   \item{U:}{ an array of dimension
#' #'   \code{TTx(DD*num_regressors_per_component)xNN}
#' #'   containing the U-type regressors; for each
#' #'   component \code{DD}; there might be a different
#' #'   number of regressors used so the number of
#' #'   columns has the component
#' #'   \code{num_regressors_per_component}; there are
#' #'   \code{NN} different cross sectional units (US
#' #'   states), so the final array dimension is
#' #'   \code{NN}}
#' #'   \item{NN:}{number of cross sectional
#' #'   units (US states)}
#' #'   \item{TT:}{number of observed time periods (e.g.
#' #'   1960-2014)}
#' #'   \item{DD:}{number of energy
#' #'   components (corresponding to dirichlet shares or
#' #'   raw counts of energy as e.g. with the
#' #'   dirichlet-multinomial)} }
#' get_data = function(data_set,
#'                     Y_to_use,
#'                     Z_to_use,
#'                     U_to_use) {
#'
#' },
