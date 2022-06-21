
./matlab_code
includes the MATLAB files for generating the results presented in the manuscript.

. MainFigure1.m: script for generating the results shown in Fig. 1.

. MainFigure2.m: script for generating the results shown in Fig. 2.

. MainFigure3.m: script for generating the results shown in Fig. 3.

. /res_opt: folder containing the optimized results for generating Figs. 1-3, to be included in /matlab_code

. covid_model_SIHM_distu.m: function used for defining the ode system implementing the SIHM model with beta_f=u*(1+w)

. covid_model_SIHM_distu_opt.m: function used for defining the ode system implementing the SIHM model with beta_f=u*(1+w)*(1-Delta_u)

. my_cost_fun_covid19_model_SIHM_distu_country: cost function for fitting procedure

. my_cost_fun_covid19_model_SIHM_distu_country_plot: function used for plotting the behavior of SIHM model for a given parameter set 

. my_cost_fun_covid19_model_SIHM_distu_opt_ph3: cost function used for implementing a control strategy with the aim to keep the peak value below the desired threshold

. data_dH_Germany.m: script for loading Germany data 

. data_dH_France.m: script for loading France data

. data_dH_Italy.m: script for loading Italy data

. data_dH_UK.m: script for loading UK data
