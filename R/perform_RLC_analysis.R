#' Perform RLC Analysis
#'
#' @param data1_l Takes a data.frame with columns id, PAR, and rETR. id is the unique RLC ID. PAR is the saturating light intensity. rETR is the relative electron transport rate
#' @param df_x A data frame used to predict y-values using the range of light intensities in column PAR e.g data.frame(PAR = seq(0.1, 926, length = 100))
#'
#' @return A list containing p0, a plot of all the fits, and final_df, containing the extracted parameters such as Ek.
#' @export
#'

#source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1")  #set theme in code
#source("https://raw.githubusercontent.com/gerard-ricardo/data/master/ssplattmy")  #for starting values

#' @export
perform_RLC_analysis <- function(data1_l, df_x) {
  # Handle small values in PAR and rETR
  data1_l$PAR <- ifelse(data1_l$PAR <= 0, 0.1, data1_l$PAR)
  data1_l$rETR <- ifelse(data1_l$rETR <= 0, 0.01, data1_l$rETR)

  # Get initial values for parameters
  starts <- data1_l %>%
    group_by(id) %>%
    do(broom::tidy(stats::getInitial(rETR ~ SSPlatt.my(PAR, alpha, beta, Ys), data = .))) %>%
    pivot_wider(names_from = names, values_from = x, names_prefix = "") %>%
    dplyr::select(., -c('NA'))
  colnames(starts) <- c("id", "alpha.s", 'beta.s', 'Pmax.s')
  starts <- NaRV.omit(starts)  # Remove infinities

  # Fit models
  fits <- data1_l  %>%
    right_join(., starts, by = 'id') %>%
    group_by(id) %>%
    do(model = try(nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)),
                         start = list(Pmax = mean(.$Pmax.s),
                                      alpha = mean(.$alpha.s),
                                      beta = mean(.$beta.s)),
                         data = .), silent = TRUE))

  # Generate predictions for each model
  usq <- list()
  for (i in 1:nrow(fits)) {
    out <- try(predict(fits$model[[i]], df_x))
    usq <- c(usq, list(out))
  }

  # Convert predictions to a data.frame
  df3 <- data.frame(t(do.call(rbind.data.frame, usq)), row.names = paste0("", 1:100))
  colnames(df3) <- as.factor(as.character(fits$id))
  df3[] <- lapply(df3, function(x) as.numeric(as.character(x)))
  df3$PAR <- df_x$PAR

  df3_long = df3 %>% pivot_longer(-PAR,  names_to = "id" ,values_to = "rETR") %>% data.frame() #keep vec.x, add all other columns to factors , add all their values to meas)
  df3_long = dplyr::arrange(df3_long, id)
  str(df3_long)
  df3_long$rETR <- as.numeric(as.character(df3_long$rETR))  #add col names
  p0 = ggplot() + geom_point(data1_l, mapping = aes(x = PAR, y = rETR), size = 1 )
  p0 = p0 + geom_line(df3_long, mapping = aes(x = PAR, y = rETR))
  p0 = p0 + facet_wrap(~id)
  p0  #clean multi-plot

  ####extract  parameters from all RLC######
  params = data1_l  %>% right_join(., starts, by = 'id') %>%
    group_by(id) %>%
    do(model = try(broom::tidy(nlsLM(rETR ~ Pmax*(1-exp(-alpha*PAR/Pmax))*(exp(-beta*PAR/Pmax)),
                                     start = list(Pmax = mean(.$Pmax.s),
                                                  alpha = mean(.$alpha.s),
                                                  beta = mean(.$beta.s)),
                                     data = .),silent = T)) )  #this get parameters for all models

  params$model[[1]]  #check for model 1
  params$model

  params$len <- sapply(params$model, length)  #using length to find non erroneous fits
  fits  = params[params$len > 2,]
  errors = params[params$len < 2,]
  unest.params = fits %>% tidyr::unnest(., model)

  df.param  = dplyr::select(unest.params, c(id, term, estimate))
  dat_wide <- df.param %>% pivot_wider(names_from = term, values_from = estimate)  #%>% dplyr::select(.,-c("NA")) #year goes to columns, their areas go as the values, area is the prefix
  dat_wide$ETRm = dat_wide$Pmax*(dat_wide$alpha/(dat_wide$alpha+dat_wide$beta))*((dat_wide$beta/(dat_wide$alpha+dat_wide$beta)))^(dat_wide$beta/dat_wide$alpha)
  dat_wide$Ek = dat_wide$ETRm/dat_wide$alpha
  dat_wide$Em =(dat_wide$Pmax/dat_wide$alpha) * log((dat_wide$alpha+dat_wide$beta)/dat_wide$beta)
  final_df = left_join(dat_wide, data1_long, by = 'id')
  final_df

  # Prepare output
  output <- list(p0 = p0, final_df = final_df)
  return(output)
}
