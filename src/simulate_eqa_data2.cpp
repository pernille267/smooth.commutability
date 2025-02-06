#include <Rcpp.h>
using namespace Rcpp;

//' Simulation of EQA data based on study design and potential differences in non-selectivity
//'
//' @title Simulation of EQA data based on study design and potential differences in non-selectivity
//' @name simulate_eqa_data2
//'
//' @param parameters A \code{list} of parametedrs used to simulate the EQA data. You must at least specify one parameter for this function to run. Except that one mandatory parameter, you may optionally choose the remaining of the parameters. These are the optimal parameters that you may include into the list:
//' \itemize{
//'   \item \code{n:} The number of samples.
//'   \item \code{R:} The number of replicates on each sample.
//'   \item \code{cvx:} The repeatability coefficient of variation for IVD-MD \code{MP_B}.
//'   \item \code{cvy:} The repeatability coefficient of variation for IVD-MD \code{MP_A}.
//'   \item \code{cil:} The lower bound of the concentration interval
//'   \item \code{ciu:} The upper bound of the concentration interval
//'   \item \code{dist:} The distribution to simulate latent variables from.
//'                      Possbile choices include \code{unif} (uniform distribution, default),
//'                      \code{norm} (normal distribution), \code{lst} (location-scale t-distribution),
//'                      \code{lnorm} (log-normal distribution)
//'   \item \code{df_tau:} The degrees of freedom for the 'lst' distribution if the distribution of latent variables are location-scale t-distributed ('lst'). Defaults to 5 if not otherwise is specified.
//'   \item \code{eta:} The heteroscedasticity factor.
//'   \item \code{eta0:} The proportion of base MS standard deviations.
//'   \item \code{qpos:} Position of systematic differences in non-selectivity. 0 signify lower range and 1 signify upper range
//'   \item \code{qran:} Interquantile range where systematic differences in non-selectivity should have its effect
//'   \item \code{prop:} average proportion of clinical samples affected by random differences in non-selectivity
//'   \item \code{mmax:} The maximum relocation magnitude in number of analytical SDs of y measurements. This assumes either prop or qpos and qran to be specified as well
//'   \item \code{b0:} For systematic linear DINS between IVD-MDs. Intercept. Defaults to 0.
//'   \item \code{b1:} For systematic linear DINS between IVD-MDs. Slope. Defaults to 1.
//'   \item \code{c0:} For systematic linear non-selectivity in IVD-MD 1. Intercept. Defaults to 0.
//'   \item \code{c1:} For systematic linear non-selectivity in IVD-MD 1. Slope. Defaults to 1.
//'   \item \code{error_dist:} The distribution to simulate measurement error components from. Possible choices include 'norm' (normal distribution, default) and 'lt' (location t-distribution)
//'   \item \code{dfx:} The degrees of freedom for the measurement error components in IVD-MD 1 if error_dist = 'lt'. Defaults to 5 if not otherwise is specified.
//'   \item \code{dfy:} The degrees of freedom for the measurement error components in IVD-MD 2 if error_dist = 'lt'. Defaults to 5 if not otherwise is specified.
//'   \item \code{md_method:} Method for simulation missing data. Possible choices include 'none' (no missing data is simulated, default), 'mar' (missing at random), 'mnar' (missing not at random) and 'marmnar' (missing both at random and not at random)
//'   \item \code{mar_prob:} The probability (value between 0 and 1) of having a missing measurement. Only relevant if \code{md_method} is 'mar' or 'marmnar'. If not specified, but \code{md_method} = 'mar' or \code{md_method} = 'marmnar', it defaults to 0.05.
//'   \item \code{mnar_threshold:} The lower bound threshold (a real value) for when a measurement should be missing. Only relevant if \code{md_method} is 'mnar' or 'marmnar'. If not specified, but \code{md_method} = 'mnar' or \code{md_method} = 'marmnar', it defaults to \code{cil}. Alternatively, if not specified, but \code{md_method} = 'mnar0' or \code{md_method} = 'marmnar0', it defaults to 0.
//' }
//' @param type \code{integer}. Set to \code{0} for default simulation of data. Set to \code{1}, \code{2} or \code{3} to simulate from custom built in non-linear functions.
//' @param AR \code{logical}. If \code{TRUE}, data is simulated including replicated measurements. Otherwise, mean of replicated measurements are returned (MOR).
//' @param include_parameters \code{logical}. If \code{TRUE}, the used parameters in the data simulation is saved and placed in a seperate list as part of the output.
//' @param shift \code{logical}. If \code{TRUE}, the simulated data change roles. MP_A becomes MP_B, and MP_B becomes MP_A.
//'
//' @description Simulates a data set with n x R rows, and four columns. The two first columns are the base ID columns (\code{SampleID} and \code{ReplicateID}). The remaining columns are numeric columns holding measurement results from the two IVD-MDs in comparison (denoted 'MP_A' (y) and 'MP_B' (x)).
//'              The form of the simulated data depends greatly on which parameters are specified in the the \code{parameters} argument.
//' @details Simulates method comparison data for \code{n} samples (e.g., clinical samples, pools, external quality assessment samples, reference material samples), where each sample is measured \code{R} times (replicated measurements). In theory, we simulate from (x_ir, y_ir) where x_ir = f(tau_i) + h_ir and y_ir = g(f(tau_i)) + v_ir.
//'          The form of f is specified through parameters \code{c0} and \code{c1}, whereas g is specified through numerous parameters such as \code{b0}, \code{b1}, \code{qpos}, \code{qran}, \code{mmax}, \code{prop}. tau_i is modelled through \code{cil}, \code{ciu}, \code{dist} and \code{df_tau}.
//'          h_ir and v_ir are measurement error components modelled through \code{cvx}, \code{cvy}. \code{cvx}, \code{cvy} can be functions of \code{error_dist}, \code{dfx}, \code{dfy}, \code{eta} and \code{eta0}.
//'          In order to convert the outputted list to a table, use either \code{as.data.frame()}, \code{as.data.table()}, \code{as.tibble()}. The most efficient way to convert is \code{setDT()} from the \code{data.table} package.
//'
//' @return A list where each list element is a column of the generated EQA data
//'
//' @examples \dontrun{
//'
//'   # Load data.table package from library
//'   library(data.table)
//'
//'   # Simulate 25 clinical samples measured in triplicate affected by
//'   # random differences in non-selectivity
//'   parameters_css <- list(n = 25, R = 3, prop = 0.1, mmax = 5, cil = 25, ciu = 75)
//'   simulated_css <- simulate_eqa_data(parameters = parameters_css)
//'
//'   # Simulate 3 external quality assessment material samples
//'   # measured in duplicate not affected by differences in non-selectivity
//'   parameters_eqams <- list(n = 3, R = 2, b0 = 0.1, b1 = 1.1)
//'   simulated_eqams <- simulate_eqa_data(parameters = parameters_eqams)
//'
//'   # We can assume that tau_i ~ lst(df_tau = 10, mu_tau = 50, var_tau = 78.583)
//'   parameters_css <- c(parameters_css, dist = "lst", df_tau = 10)
//'   simulated_css_lst <- simulate_eqa_data(parameters = parameters_css)
//'
//'   # We can convert the list objects to data.table objects using setDT()
//'   setDT(simulated_css)
//'   setDT(simulated_eqams)
//'   setDT(simulated_css_lst)
//'
//'   # Print results
//'   print(simulated_css)
//'   print(simulated_eqams)
//'   print(simulated_css_lst)
//'
//' }
//'


// [[Rcpp::export]]
List simulate_eqa_data2(List parameters, int type = 1, bool AR = false, bool include_parameters = false, bool shift = false){

 // Assume that no parameters are given as default
 int n_exists = 0;
 int R_exists = 0;
 int cvx_exists = 0;
 int cvy_exists = 0;
 int cil_exists = 0;
 int ciu_exists = 0;
 int dist_exists = 0;
 int df_tau_exists = 0;
 int eta_exists = 0;
 int eta0_exists = 0;
 int qpos_exists = 0;
 int qran_exists = 0;
 int prop_exists = 0;
 int mmax_exists = 0;
 int b0_exists = 0;
 int b1_exists = 0;
 int c0_exists = 0;
 int c1_exists = 0;
 int error_dist_exists = 0;
 int dfx_exists = 0;
 int dfy_exists = 0;
 int md_method_exists = 0;
 int mar_prob_exists = 0;
 int mnar_threshold_exists = 0;
 int qdir_exists = 0;
 int obs_tau_exists = 0;
 int qlim_exists = 0;
 int qnor_exists = 0;
 int sdx_exists = 0;
 int sdy_exists = 0;

 // Checks which of the parameters found in the 'parameters' argument
 CharacterVector given_parameters = parameters.names();
 int number_given_parameters = given_parameters.size();

 for(int i = 0; i < number_given_parameters; ++i){

   CharacterVector candidate_param(30);
   candidate_param[0] = "n";
   candidate_param[1] = "R";
   candidate_param[2] = "cvx";
   candidate_param[3] = "cvy";
   candidate_param[4] = "cil";
   candidate_param[5] = "ciu";
   candidate_param[6] = "dist";
   candidate_param[7] = "df_tau";
   candidate_param[8] = "eta";
   candidate_param[9] = "eta0";
   candidate_param[10] = "qpos";
   candidate_param[11] = "qran";
   candidate_param[12] = "prop";
   candidate_param[13] = "mmax";
   candidate_param[14] = "b0";
   candidate_param[15] = "b1";
   candidate_param[16] = "c0";
   candidate_param[17] = "c1";
   candidate_param[18] = "error_dist";
   candidate_param[19] = "dfx";
   candidate_param[20] = "dfy";
   candidate_param[21] = "md_method";
   candidate_param[22] = "mar_prob";
   candidate_param[23] = "mnar_threshold";
   candidate_param[24] = "qdir";
   candidate_param[25] = "obs_tau";
   candidate_param[26] = "qlim";
   candidate_param[27] = "qnor";
   candidate_param[28] = "sdx";
   candidate_param[29] = "sdy";



   if(candidate_param[0] == given_parameters[i]){
     n_exists = 1;
   }
   else if(candidate_param[1] == given_parameters[i]){
     R_exists = 1;
   }
   else if(candidate_param[2] == given_parameters[i]){
     cvx_exists = 1;
   }
   else if(candidate_param[3] == given_parameters[i]){
     cvy_exists = 1;
   }
   else if(candidate_param[4] == given_parameters[i]){
     cil_exists = 1;
   }
   else if(candidate_param[5] == given_parameters[i]){
     ciu_exists = 1;
   }
   else if(candidate_param[6] == given_parameters[i]){
     dist_exists = 1;
   }
   else if(candidate_param[7] == given_parameters[i]){
     df_tau_exists = 1;
   }
   else if(candidate_param[8] == given_parameters[i]){
     eta_exists = 1;
   }
   else if(candidate_param[9] == given_parameters[i]){
     eta0_exists = 1;
   }
   else if(candidate_param[10] == given_parameters[i]){
     qpos_exists = 1;
   }
   else if(candidate_param[11] == given_parameters[i]){
     qran_exists = 1;
   }
   else if(candidate_param[12] == given_parameters[i]){
     prop_exists = 1;
   }
   else if(candidate_param[13] == given_parameters[i]){
     mmax_exists = 1;
   }
   else if(candidate_param[14] == given_parameters[i]){
     b0_exists = 1;
   }
   else if(candidate_param[15] == given_parameters[i]){
     b1_exists = 1;
   }
   else if(candidate_param[16] == given_parameters[i]){
     c0_exists = 1;
   }
   else if(candidate_param[17] == given_parameters[i]){
     c1_exists = 1;
   }
   else if(candidate_param[18] == given_parameters[i]){
     error_dist_exists = 1;
   }
   else if(candidate_param[19] == given_parameters[i]){
     dfx_exists = 1;
   }
   else if(candidate_param[20] == given_parameters[i]){
     dfy_exists = 1;
   }
   else if(candidate_param[21] == given_parameters[i]){
     md_method_exists = 1;
   }
   else if(candidate_param[22] == given_parameters[i]){
     mar_prob_exists = 1;
   }
   else if(candidate_param[23] == given_parameters[i]){
     mnar_threshold_exists = 1;
   }
   else if(candidate_param[24] == given_parameters[i]){
     qdir_exists = 1;
   }
   else if(candidate_param[25] == given_parameters[i]){
     obs_tau_exists = 1;
   }
   else if(candidate_param[26] == given_parameters[i]){
     qlim_exists = 1;
   }
   else if(candidate_param[27] == given_parameters[i]){
     qnor_exists = 1;
   }
   else if(candidate_param[28] == given_parameters[i]){
     sdx_exists = 1;
   }
   else if(candidate_param[29] == given_parameters[i]){
     sdy_exists = 1;
   }
 }

 // Base parameters
 int n = 25;
 int R = 3;
 double cvx = 0;
 double cvy = 0;
 double cil = 0;
 double ciu = 0;

 // Heteroscedasticity parameters
 double eta = 0;
 double eta0 = 0;

 // Random and systematic differences in non-selectivity parameters
 int qpos = 0;
 int qdir = 0;
 double qran = 0;
 double qlim = 0;
 double qnor = 1;
 double prop = 0;
 double mmax = 0;

 // Systematic linear differences in non-selectivity parameters
 double b0 = 0;
 double b1 = 1;

 // Systematic linear non-selectivity
 double c0 = 0;
 double c1 = 1;

 // Modelling assumptions
 std::string dist = "unif";
 std::string error_dist = "norm";

 // Parameters for modelling assumptions
 double sigma_tau = 0;
 double mu_tau = 0;
 double df_tau = 5.0;
 double dfx = 5.0;
 double dfy = 5.0;

 // Parameters for modelling missing data
 std::string md_method = "none";
 double mar_prob = 0;
 double mnar_threshold = 0;

 // Generation of n
 if(n_exists == 1){
   int reg_n = parameters["n"];
   n = reg_n;
 }
 else{
   int gen_n = R::rpois(25);
   IntegerVector ns (2);
   ns[0] = gen_n;
   ns[1] = 30;
   gen_n = min(ns);
   ns[0] = gen_n;
   ns[1] = 20;
   n = max(ns);
 }
 // Generation of R
 if(R_exists == 1){
   int reg_R = parameters["R"];
   R = reg_R;
 }
 else{
   double u = R::runif(0,1);
   if(u < 0.10){
     R = 2;
   }
   else if(u < 0.95){
     R = 3;
   }
   else{
     R = 4;
   }
 }

 if(cil_exists == 1){
   double reg_cil = parameters["cil"];
   cil = reg_cil;
 }
 else{
   if(type != 0){
     cil = 2;
   }
   else{
     cil = R::rf(1.057057, 8.15) * 44;
     if(cil < 0.01){
       cil = R::rf(1.057057, 8.15) * 44;
       if(cil < 0.01){
         cil = R::rf(1.057057, 8.15) * 44;
         if(cil < 0.01){
           cil = R::rf(1.057057, 8.15) * 44;
         }
       }
     }
   }
 }
 if(ciu_exists == 1){
   double reg_ciu = parameters["ciu"];
   if(reg_ciu <= cil){
     Rcout << "ciu is" << reg_ciu << "which is smaller than" << cil << "\n";
     stop("ciu must be larger than cil!");
   }
   ciu = reg_ciu;
 }
 else{
   if(cil > 0){
     double multiplier = R::rbeta(0.78, 11) * 44;
     if(type != 0){
       ciu += 10.0;
     }
     else{
       ciu += cil + cil * multiplier;
     }

   }
   else{
     stop("cil is negative or zero. This should not be possible! Check if something is wrong");
   }
 }

 // Generation of cvx
 if(cvx_exists == 1){
   double reg_cvx = parameters["cvx"];
   cvx += reg_cvx;
 }
 else{
   if(sdx_exists == 1){
     double reg_sdx = parameters["sdx"];
     cvx += reg_sdx / (0.5 * (cil + ciu));
   }
   else{
     cvx += R::rbeta(2, 5) / 10.0;
   }
 }
 // Generation of cvy
 if(cvy_exists == 1){
   double reg_cvy = parameters["cvy"];
   cvy += reg_cvy;
 }
 else{
   if(sdy_exists == 1){
     double reg_sdy = parameters["sdy"];
     cvy += reg_sdy / (0.5 * (cil + ciu));
   }
   else{
     cvy += R::rbeta(2, 5) / 10.0;
   }
 }

 if(dist_exists == 1){
   std::string reg_dist = parameters["dist"];
   dist = reg_dist;
   if(dist == "norm"){
     mu_tau += 0.5 * (cil + ciu);
     sigma_tau += 0.2149292 * (ciu - cil);
   }
   else if(dist == "lst"){
     if(df_tau_exists == 1){
       double reg_df_tau = parameters["df_tau"];
       df_tau = reg_df_tau;
     }
     mu_tau += 0.5 * (cil + ciu);
     sigma_tau += 0.49 * (ciu - cil) * (1.0 / R::qt(0.99, df_tau, 0, 0));
   }
   else if(dist == "lnorm"){
     mu_tau += 0.5 * (log(cil + 0.99 * (ciu - cil)) + log(cil + 0.01 * (ciu - cil)));
     sigma_tau += 0.2149292 * log((cil + 0.99 * (ciu - cil))/(cil + 0.01 * (ciu - cil)));
   }
   else{
     mu_tau += 0.5 * (cil + ciu);
     sigma_tau += (1.0 / sqrt(12.0)) * (ciu - cil);
   }
 }

 if(error_dist_exists == 1){
   std::string reg_error_dist = parameters["error_dist"];
   error_dist = reg_error_dist;
   if(error_dist != "norm" and (error_dist == "lt" or error_dist == "t" or error_dist == "lst")){
     error_dist = "lt";
     if(dfx_exists == 1){
       double reg_dfx = parameters["dfx"];
       dfx = reg_dfx;
     }
     if(dfy_exists == 1){
       double reg_dfy = parameters["dfy"];
       dfy = reg_dfy;
     }
   }
 }

 if(eta_exists == 1){
   double reg_eta = parameters["eta"];
   eta += reg_eta;
 }
 else{
   eta = 1;
 }
 if(eta0_exists == 1){
   double reg_eta0 = parameters["eta0"];
   eta0 = reg_eta0;
 }
 else{
   eta0 = 1;
 }

 if(b0_exists == 1){
   double reg_b0 = parameters["b0"];
   b0 = reg_b0;
 }

 if(b1_exists == 1){
   double reg_b1 = parameters["b1"];
   b1 = reg_b1;
 }

 if(c0_exists == 1){
   double reg_c0 = parameters["c0"];
   c0 = reg_c0;
 }

 if(c1_exists == 1){
   double reg_c1 = parameters["c1"];
   c1 = reg_c1;
 }

 if(qran_exists == 1 and prop_exists == 1){
   //Rcout << "Both qran and prop are found in parameters" << "\n";
   //Rcout << "Make sure only one of them is used" << "\n";
   //Rcout << "prop is removed in favor of qran" << "\n";
   prop_exists -= 1;
 }
 if(qpos_exists == 1){
   double reg_qpos = parameters["qpos"];
   qpos += reg_qpos;
 }
 else{
   qpos -= 1;
 }
 if(qran_exists == 1){
   double reg_qran = parameters["qran"];
   qran += reg_qran;
 }
 else{
   qran += 0;
 }

 if(qdir_exists == 1){
   int reg_qdir = parameters["qdir"];
   qdir = reg_qdir;
   if(qdir >= 1){
     qdir = 1;
   }
   else if(qdir <= -1){
     qdir = -1;
   }
   else{
     qdir = 1;
   }
 }
 else{
   int above = R::rbinom(1, 0.5);
   if(above == 1){
     qdir += 1;
   }
   else{
     qdir -= 1;
   }
 }

 if(prop_exists == 1){
   double reg_prop = parameters["prop"];
   prop += reg_prop;
 }
 else{
   prop = 0;
 }
 if(mmax_exists == 1){
   double reg_mmax = parameters["mmax"];
   mmax += reg_mmax;
 }

 if(md_method_exists == 1){
   std::string reg_md_method = parameters["md_method"];
   md_method = reg_md_method;

   if(md_method == "mar"){
     if(mar_prob_exists == 1){
       double reg_mar_prob = parameters["mar_prob"];
       mar_prob += reg_mar_prob;
     }
     else{
       mar_prob += 0.05;
     }
   }
   else if(md_method == "mnar"){
     if(mnar_threshold_exists == 1){
       double reg_mnar_threshold = parameters["mnar_threshold"];
       mnar_threshold += reg_mnar_threshold;
     }
     else{
       mnar_threshold += cil;
     }
   }
   else if(md_method == "mnar0"){
     if(mnar_threshold_exists == 1){
       double reg_mnar_threshold = parameters["mnar_threshold"];
       mnar_threshold += reg_mnar_threshold;
     }
   }
   else if(md_method == "marmnar"){
     if(mnar_threshold_exists == 1){
       double reg_mnar_threshold = parameters["mnar_threshold"];
       mnar_threshold += reg_mnar_threshold;
     }
     else{
       mnar_threshold += cil;
     }
     if(mar_prob_exists == 1){
       double reg_mar_prob = parameters["mar_prob"];
       mar_prob += reg_mar_prob;
     }
     else{
       mar_prob += 0.05;
     }
   }

   else if(md_method == "marmnar0"){
     if(mnar_threshold_exists == 1){
       double reg_mnar_threshold = parameters["mnar_threshold"];
       mnar_threshold += reg_mnar_threshold;
     }
     if(mar_prob_exists == 1){
       double reg_mar_prob = parameters["mar_prob"];
       mar_prob += reg_mar_prob;
     }
     else{
       mar_prob += 0.05;
     }
   }
   else{
     mnar_threshold -= 999999.99;
   }
 }
 else{
   mnar_threshold -= 999999.99;
 }

 if(obs_tau_exists == 1){
   NumericVector reg_obs_tau = parameters["obs_tau"];
   n = reg_obs_tau.size();
   dist = "none";
 }

 NumericVector unsorted_tau(n);
 NumericVector tau(n);

 int nR = n * R;
 IntegerVector SampleID(nR);
 IntegerVector SampleID2(n);
 IntegerVector ReplicateID(nR);
 NumericVector MP_A(nR);
 NumericVector MP_A2(n);
 NumericVector MP_B(nR);
 NumericVector MP_B2(n);

 double average = 0.5 * (cil + ciu);
 if(dist == "lnorm"){
   average = exp(mu_tau + pow(sigma_tau, 2.0) / 2.0);
 }

 double base_x = average * cvx;
 double base_y = average * cvy;

 double beg_sdx = base_x * eta0;
 double end_sdx = base_x * eta * eta0;
 double seg_sdx = 0;

 if(eta0 > 0 and eta > 0){
   seg_sdx = (end_sdx - beg_sdx) / n;
 }

 double beg_sdy = base_y * eta0;
 double end_sdy = base_y * eta * eta0;
 double seg_sdy = 0;

 if(eta0 > 0 and eta > 0){
   seg_sdy = (end_sdy - beg_sdy) / n;
 }

 if(dist == "none"){
   NumericVector reg_obs_tau = parameters["obs_tau"];
   for(int i = 0; i < n; ++i){
     unsorted_tau[i] = reg_obs_tau[i];
   }
 }
 else if(dist == "norm"){
   for(int i = 0; i < n; ++i){
     unsorted_tau[i] = R::rnorm(mu_tau, sigma_tau);
   }
 }
 else if(dist == "lst"){
   for(int i = 0; i < n; ++i){
     unsorted_tau[i] = mu_tau + sigma_tau * R::rt(df_tau);
   }
 }
 else if(dist == "lnorm"){
   for(int i = 0; i < n; ++i){
     unsorted_tau[i] = R::rlnorm(mu_tau, sigma_tau);
   }
 }
 else{
   for(int i = 0; i < n; ++i){
     unsorted_tau[i] = R::runif(cil, ciu);
   }
 }


 if(eta0_exists == 1 and eta_exists == 1){
   tau = unsorted_tau.sort();
 }

 for(int i = 0; i < n; ++i){
   SampleID2[i] = i + 1;
   int na_count_A = 0;
   int na_count_B = 0;
   double hetero_extra_x = i * seg_sdx;
   double hetero_extra_y = i * seg_sdy;
   double sdx = beg_sdx + hetero_extra_x;
   double sdy = beg_sdy + hetero_extra_y;

   if(sdx <= 0 or sdy <= 0){
     stop("Your choices of eta and eta0 resulted in negative standard deviations. Calculations are terminated");
   }

   if(eta0_exists == 0 or eta_exists == 0){
     tau[i] = unsorted_tau[i];
   }

   double limit = 0;
   double relocating_magnitude = 0;
   int relocate_sample_i = 0;

   // if prop exists, and CS is dins-affected, randomly select relocation magnitude from from Beta(2, 2) * mmax
   if(prop_exists == 1){
     relocate_sample_i += R::rbinom(1, prop);
     double sign_relocated_sample_i = R::runif(0, 1);
     if(sign_relocated_sample_i < 0.5){
       sign_relocated_sample_i = -1;
     }
     else{
       sign_relocated_sample_i = 1;
     }
     relocating_magnitude += R::rbeta(2, 2) * sign_relocated_sample_i * mmax;
   }
   else if(qran_exists == 1){
     if(qpos == 0){
       if(qlim_exists == 1){
         double reg_qlim = parameters["qlim"];
         limit += reg_qlim;
       }
       else{
         limit += cil + qran * (ciu - cil);
       }
       if(qnor_exists == 1){
         double reg_qnor = parameters["qnor"];
         qnor = reg_qnor;
       }
       else{
         qnor = limit - cil;
       }
       qlim = limit;
       if(tau[i] <= limit){
         relocate_sample_i = relocate_sample_i + 1;
         //double rel_diff = (limit - tau[i]) / (limit - cil);
         double rel_diff = R::pbeta((limit - tau[i]) / qnor, 2, 2, 1, 0) / 2.0;
         //relocating_magnitude = rel_diff * mmax * qdir;
         relocating_magnitude = 2 * rel_diff * mmax * qdir;
       }
     }
     else if(qpos == 1){
       if(qlim_exists == 1){
         double reg_qlim = parameters["qlim"];
         limit += reg_qlim;
       }
       else{
         limit += ciu - qran * (ciu - cil);
       }
       if(qnor_exists == 1){
         double reg_qnor = parameters["qnor"];
         qnor = reg_qnor;
       }
       else{
         qnor = ciu - limit;
       }
       qlim = limit;
       if(tau[i] >= limit){
         relocate_sample_i = relocate_sample_i + 1;
         //double rel_diff = (tau[i] - limit) / (ciu - limit);
         double rel_diff = R::pbeta((tau[i] - limit) / qnor, 2, 2, 1, 0) / 2.0;
         //relocating_magnitude = rel_diff * mmax * qdir;
         relocating_magnitude = 2 * rel_diff * mmax * qdir;
       }
     }
   }
   int lower = i * R;
   int upper = R * (1 + i) - 1;
   IntegerVector idx = Rcpp::Range(lower, upper);
   if(error_dist == "norm"){
     for(int r = 0; r < R; ++r){
       SampleID[idx[r]] = i + 1;
       ReplicateID[idx[r]] = r + 1;
       MP_B[idx[r]] = c0 + tau[i] * c1;
       if(type == 1){
         MP_A[idx[r]] = MP_B[idx[r]] + 0.9 * sin(0.4 * pow(MP_B[idx[r]], 1.06)) + R::rnorm(0, sdy);
       }
       else if(type == 2){
         MP_A[idx[r]] = MP_B[idx[r]] + 0.05 * exp(0.16 * pow(MP_B[idx[r]], 1.35)) + R::rnorm(0, sdy);
       }
       else if(type == 3){
         MP_A[idx[r]] = MP_B[idx[r]] - exp(-0.5 * pow((MP_B[idx[r]] - 1.5)/2.0, 2)) + R::rnorm(0, sdy);
       }
       else{
         MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + R::rnorm(0, sdy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
       }
       MP_B[idx[r]] = MP_B[idx[r]] + R::rnorm(0, sdx);
       if(!AR){
         MP_A2[i] += MP_A[idx[r]];
         MP_B2[i] += MP_B[idx[r]];
       }
       if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
         bool mar_r_A = R::runif(0, 1) < mar_prob;
         bool mar_r_B = R::runif(0, 1) < mar_prob;
         bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
         bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
         if(mar_r_A or mnar_r_A){
           MP_A[idx[r]] = NA_REAL;
           if(!AR){
             na_count_A += 1;
             MP_A2[i] -= MP_A[idx[r]];
           }

         }
         if(mar_r_B or mnar_r_B){
           MP_B[idx[r]] = NA_REAL;
           if(!AR){
             na_count_B += 1;
             MP_B2[i] -= MP_B[idx[r]];
           }
         }
       }
       else{
         if(MP_A[idx[r]] < 0){
           MP_A[idx[r]] = MP_A[idx[r]] * (-1);
         }
         if(MP_B[idx[r]] < 0){
           MP_B[idx[r]] = MP_B[idx[r]] * (-1);
         }
       }
     }
   }
   else if(error_dist == "lt"){
     for(int r = 0; r < R; ++r){
       SampleID[idx[r]] = i + 1;
       ReplicateID[idx[r]] = r + 1;
       MP_B[idx[r]] = c0 + tau[i] * c1;
       if(type == 1){
         MP_A[idx[r]] = MP_B[idx[r]] + 0.9 * sin(0.4 * pow(MP_B[idx[r]], 1.06)) + sdy * R::rt(dfy);
       }
       else if(type == 2){
         MP_A[idx[r]] = MP_B[idx[r]] + 0.05 * exp(0.16 * pow(MP_B[idx[r]], 1.35)) + sdy * R::rt(dfy);
       }
       else if(type == 3){
         MP_A[idx[r]] = MP_B[idx[r]] - exp(-0.5 * pow((MP_B[idx[r]] - 1.5)/2.0, 2)) + sdy * R::rt(dfy);
       }
       else{
         MP_A[idx[r]] = b0 + MP_B[idx[r]] * b1 + sdy * R::rt(dfy) + relocate_sample_i * relocating_magnitude * sqrt(pow(sdx, 2) + pow(sdy, 2));
       }
       MP_B[idx[r]] = MP_B[idx[r]] + sdx * R::rt(dfx);
       if(!AR){
         MP_A2[i] += MP_A[idx[r]];
         MP_B2[i] += MP_B[idx[r]];
       }
       if(md_method == "mar" or md_method == "mnar" or md_method == "mnar0" or md_method == "marmnar" or md_method == "marmnar0"){
         bool mar_r_A = R::runif(0, 1) < mar_prob;
         bool mar_r_B = R::runif(0, 1) < mar_prob;
         bool mnar_r_A = MP_A[idx[r]] < mnar_threshold;
         bool mnar_r_B = MP_B[idx[r]] < mnar_threshold;
         if(mar_r_A or mnar_r_A){
           MP_A[idx[r]] = NA_REAL;
         }
         if(mar_r_B or mnar_r_B){
           MP_B[idx[r]] = NA_REAL;
         }
       }
       else{
         if(MP_A[idx[r]] < 0){
           MP_A[idx[r]] = MP_A[idx[r]] * (-1);
         }
         if(MP_B[idx[r]] < 0){
           MP_B[idx[r]] = MP_B[idx[r]] * (-1);
         }
       }
     }
   }
   if(!AR){
     if(R - na_count_A > 0){
       MP_A2[i] = MP_A2[i] / (R - na_count_A);
     }
     else{
       MP_A2[i] = NA_REAL;
     }
     if(R - na_count_A > 0){
       MP_B2[i] = MP_B2[i] / (R - na_count_B);
     }
     else{
       MP_B2[i] = NA_REAL;
     }
   }
 }

 List used_parameters = List::create(Named("n") = n,
                                     Named("R") = R,
                                     Named("cvx") = cvx,
                                     Named("cvy") = cvy,
                                     Named("sdx") = base_x,
                                     Named("sdy") = base_y,
                                     Named("cil") = cil,
                                     Named("ciu") = ciu,
                                     Named("dfx") = dfx,
                                     Named("dfy") = dfy,
                                     Named("qpos") = qpos,
                                     Named("qran") = qran,
                                     Named("qdir") = qdir,
                                     Named("qlim") = qlim,
                                     Named("qnor") = qnor,
                                     Named("mmax") = mmax,
                                     Named("md_method") = md_method,
                                     Named("mar_prob") = mar_prob,
                                     Named("mnar_threshold") = mnar_threshold);

 if(!AR){
   if(include_parameters){
     if(shift){
       List sim_data = List::create(Named("SampleID") = SampleID2, Named("MP_A") = round(MP_B2, 6), Named("MP_B") = round(MP_A2, 6));
       List out = List::create(Named("simulated_data") = sim_data, Named("parameters") = used_parameters);
       return out;
     }
     List sim_data = List::create(Named("SampleID") = SampleID2, Named("MP_A") = round(MP_A2, 6), Named("MP_B") = round(MP_B2, 6));
     List out = List::create(Named("simulated_data") = sim_data, Named("parameters") = used_parameters);
     return out;
   }
   if(shift){
     List out = List::create(Named("SampleID") = SampleID2, Named("MP_A") = round(MP_B2, 6), Named("MP_B") = round(MP_A2, 6));
     return out;
   }
   List out = List::create(Named("SampleID") = SampleID2, Named("MP_A") = round(MP_A2, 6), Named("MP_B") = round(MP_B2, 6));
   return out;
 }
 if(include_parameters){
   if(shift){
     List sim_data = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = round(MP_B, 6), Named("MP_B") = round(MP_A, 6));
     List out = List::create(Named("simulated_data") = sim_data, Named("parameters") = used_parameters);
     return out;
   }
   List sim_data = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = round(MP_A, 6), Named("MP_B") = round(MP_B, 6));
   List out = List::create(Named("simulated_data") = sim_data, Named("parameters") = used_parameters);
   return out;
 }
 List out = List::create(Named("SampleID") = SampleID, Named("ReplicateID") = ReplicateID, Named("MP_A") = round(MP_A, 6), Named("MP_B") = round(MP_B, 6));
 return out;
}
