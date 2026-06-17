/*===========================================================================
  16_MODEL_IRF.DO
  Transmission block: map the model-implied spread path into output/credit
  responses via the log-linearized banking equations (Section 9, IR.1-IR.5),
  then (i) SMM-calibrate the deep transmission parameters {xi, phi, Phi_N}
  to the EMPIRICAL non-default IRF, and (ii) validate the DEFAULT path
  out-of-sample.

  Implements the CORRECTED recursions (iterating LL.5a and LL.3a directly,
  which avoids the closed-form IR.1 summation error flagged in the audit):
      net worth :  n_h  = Phi_N*n_{h-1} - B*s_h - Phi_O*gamma*s_{h-1}     (LL.5a)
      capital   :  k_h  = (1-delta)*k_{h-1} - eta*(gamma*s_{h-1}-Omega*n_{h-1}) (LL.3a)
      output    :  y_h  = alpha*k_h - eps_p*(gamma*s_h - Omega*n_h)        (LL.2a)
      credit    :  l_h  = (lambda*N_ss/l_ss)*n_h                           (IR.4)

  REQUIRES: 14_calibration.do, 15_solve_default.do, and the empirical IRFs
            (irf_nd.dta, irf_def.dta from 02/03).

  Output: fig_model_vs_data.pdf  +  cal_transmission.dta (fitted params)
===========================================================================*/

* ── Load calibration + model spread event path ─────────────────────────────
use "$clean/cal_params.dta", clear
foreach p in beta alpha delta s_ss bB_gdp N_gdp credit_gdp lambda Rstar Phi_O {
    scalar `p' = `p'[1]
}

* Model-implied spread path (decimal), horizons -2..4
preserve
    use "$clean/cal_spread_event.dta", clear
    mkmat horizon s_model, matrix(SPATH)
restore

* Empirical non-default output IRF (target), horizons -2..4 -> keep 0..4
preserve
    use "$clean/irf_nd.dta", clear
    keep if horizon >= 0
    mkmat horizon b, matrix(IRF_ND)
restore

* Empirical default-linked output IRF (validation target)
preserve
    use "$clean/irf_def.dta", clear
    keep if horizon >= 0
    mkmat horizon b, matrix(IRF_DEF)
restore

mata:
mata clear

// ─── Fixed calibration ──────────────────────────────────────────────────
beta       = st_numscalar("beta");
alpha      = st_numscalar("alpha");
delta      = st_numscalar("delta");
s_ss       = st_numscalar("s_ss");
mu         = st_numscalar("mu");
bBN        = st_numscalar("bB_gdp")/st_numscalar("N_gdp");
credit_gdp = st_numscalar("credit_gdp");
N_gdp      = st_numscalar("N_gdp");
lambda     = st_numscalar("lambda");
Phi_O      = st_numscalar("Phi_O");
gamma      = st_numscalar("gamma");

// ── Helper: robust horizon match (round avoids svmat float precision issues) ──
function pickrow(M, h) {
    idx = selectindex(round(M[.,1],1) :== h);
    return(M[idx, 2]);
}

// ─── Spread path: horizons 0..4, deviation from SS ───────────────────────
SP   = st_matrix("SPATH");
shat = J(5,1,0);
for (h=0; h<=4; h++) {
    shat[h+1] = pickrow(SP, h) - s_ss;
}
printf("\n=== MODEL SPREAD DEVIATIONS (shat, decimal) ===\n");
for (h=0;h<=4;h++) printf("  h=%g  shat=%8.5f  (%4.0f bps)\n", h, shat[h+1], shat[h+1]*10000);

// ─── Empirical IRFs ───────────────────────────────────────────────────────
ND  = st_matrix("IRF_ND");
DEF = st_matrix("IRF_DEF");
bnd  = J(5,1,.);
bdef = J(5,1,.);
for (h=0;h<=4;h++) {
    bnd[h+1]  = pickrow(ND,  h);
    bdef[h+1] = pickrow(DEF, h);
}

// ─── Empirical credit IRFs ─────────────────────────────────────────────────
// irf_nd_credit.dta and irf_def_credit.dta produced by 03_lp_resolution.do
// If missing, set to zero (graceful fallback)
stata("capture confirm file \"$clean/irf_nd_credit.dta\"");
if (st_numscalar("_rc")==0) {
    stata("preserve");
    stata("use \"$clean/irf_nd_credit.dta\", clear");
    stata("mkmat horizon b, matrix(CND)");
    stata("restore");
    CND = st_matrix("CND");
    bnd_cred = J(5,1,.);
    for (h=0;h<=4;h++) bnd_cred[h+1] = pickrow(CND, h);
} else {
    bnd_cred = J(5,1,0);
    printf("  (credit ND IRF not found — credit plot skipped)\n");
}
stata("capture confirm file \"$clean/irf_def_credit.dta\"");
if (st_numscalar("_rc")==0) {
    stata("preserve");
    stata("use \"$clean/irf_def_credit.dta\", clear");
    stata("mkmat horizon b, matrix(CDEF)");
    stata("restore");
    CDEF = st_matrix("CDEF");
    bdef_cred = J(5,1,.);
    for (h=0;h<=4;h++) bdef_cred[h+1] = pickrow(CDEF, h);
} else {
    bdef_cred = J(5,1,0);
    printf("  (credit DEF IRF not found — credit plot skipped)\n");
}

// ─── Transmission function: output + net worth ─────────────────────────────
// returns (H x 2): col1=yhat (pp), col2=nhat (log-dev, NOT in pp)
function trans_full(xi, phi, Phi_N, shat, gamma, alpha, delta, beta, bBN, Phi_O, s_ss) {
    RL_ss = 1/beta + phi*bBN;
    nexp  = (1-alpha)/alpha;
    eps_p = xi*nexp*RL_ss/(1 + xi*(RL_ss-1));
    eta   = 1/((1-alpha)*(RL_ss-(1-delta)));
    Omega = phi*bBN;
    Bsens = bBN/(1+s_ss)^2;
    H = 5;
    nw   = J(H,1,0);
    kap  = J(H,1,0);
    yhat = J(H,1,0);
    // h=0 (index 1)
    nw[1]   = -Bsens*shat[1];
    kap[1]  = 0;
    yhat[1] = alpha*kap[1] - eps_p*(gamma*shat[1] - Omega*nw[1]);
    // h=1..4
    for (h=2; h<=H; h++) {
        nw[h]   = Phi_N*nw[h-1] - Bsens*shat[h] - Phi_O*gamma*shat[h-1];
        kap[h]  = (1-delta)*kap[h-1] - eta*(gamma*shat[h-1] - Omega*nw[h-1]);
        yhat[h] = alpha*kap[h] - eps_p*(gamma*shat[h] - Omega*nw[h]);
    }
    return((yhat:*100, nw));    // col1=output pp, col2=nhat
}

// ─── SMM: minimise SSE on output only (xi range extended to 0.05) ──────────
xig = range(0.05, 0.90, 0.05);
phg = range(0.01, 0.15, 0.01);
png = range(0.30, 0.95, 0.05);

best = 1e12; bx=.; bp=.; bn=.;
for (i=1;i<=rows(xig);i++) {
  for (j=1;j<=rows(phg);j++) {
    for (k=1;k<=rows(png);k++) {
      res = trans_full(xig[i], phg[j], png[k], shat, gamma, alpha, delta, beta, bBN, Phi_O, s_ss);
      yh  = res[.,1];
      sse = sum((yh - bnd):^2);
      if (sse < best) { best=sse; bx=xig[i]; bp=phg[j]; bn=png[k]; }
    }
  }
}
printf("\n=== SMM FIT (non-default output IRF) ===\n");
printf("  xi*    = %5.3f\n", bx);
printf("  phi*   = %5.3f\n", bp);
printf("  Phi_N* = %5.3f\n", bn);
printf("  SSE    = %8.4f\n", best);

// ─── Final non-default IRFs at fitted params ───────────────────────────────
res_nd  = trans_full(bx, bp, bn, shat, gamma, alpha, delta, beta, bBN, Phi_O, s_ss);
yfit    = res_nd[.,1];
nfit    = res_nd[.,2];
// private credit: lhat = (lambda*N_ss/ell_ss)*nhat  (in pp)
lev_amp = lambda * N_gdp / credit_gdp;
lfit    = nfit :* lev_amp :* 100;

printf("\n  h   data_y(nd)  model_y   data_l(nd)  model_l\n");
for (h=0;h<=4;h++) {
    printf(" %g   %7.3f    %7.3f    %7.3f     %7.3f\n",
        h, bnd[h+1], yfit[h+1], bnd_cred[h+1], lfit[h+1]);
}

st_matrix("YFIT",  yfit);
st_matrix("LFIT",  lfit);
st_numscalar("xi_fit",   bx);
st_numscalar("phi_fit",  bp);
st_numscalar("PhiN_fit", bn);

// ─── Default path (out-of-sample) ─────────────────────────────────────────
RL_ss  = 1/beta + bp*bBN;
nexp   = (1-alpha)/alpha;
eps_p  = bx*nexp*RL_ss/(1+bx*(RL_ss-1));
dRL_aut = -bdef[1]/(eps_p*100);
printf("\n=== DEFAULT PATH VALIDATION ===\n");
printf("  dRL_aut (pinned to h=0) = %6.4f\n", dRL_aut);

ydef = J(5,1,.);
ldef = J(5,1,0);   // credit ≈ 0 (gambling for resurrection: ell=0, GDP normalises)
printf("\n  h   data_y(def)  model_y   [h>=1=oos]\n");
for (h=0;h<=4;h++) {
    surv        = (1-mu)^h;
    kdef        = surv * h * ln(1-delta);
    ydef[h+1]   = (alpha*kdef - eps_p*dRL_aut*surv)*100;
    tag = (h==0 ? "(pinned)" : "(oos)");
    printf(" %g   %7.3f      %7.3f   %s\n", h, bdef[h+1], ydef[h+1], tag);
}
st_matrix("YDEF", ydef);
st_matrix("LDEF", ldef);

end

* ── Assemble model-vs-data dataset ────────────────────────────────────────
clear
set obs 5
gen horizon = _n - 1

svmat YFIT, names(ymod_nd)
svmat YDEF, names(ymod_def)
svmat LFIT, names(lmod_nd)
svmat LDEF, names(lmod_def)
rename ymod_nd1  ymod_nd
rename ymod_def1 ymod_def
rename lmod_nd1  lmod_nd
rename lmod_def1 lmod_def

* merge empirical output IRFs
merge 1:1 horizon using "$clean/irf_nd.dta",  keepusing(b) nogen keep(master match)
rename b bdata_nd
merge 1:1 horizon using "$clean/irf_def.dta", keepusing(b) nogen keep(master match)
rename b bdata_def

* merge empirical credit IRFs (optional — graceful if missing)
capture merge 1:1 horizon using "$clean/irf_nd_credit.dta",  keepusing(b) nogen keep(master match)
capture rename b bdata_nd_cred
capture merge 1:1 horizon using "$clean/irf_def_credit.dta", keepusing(b) nogen keep(master match)
capture rename b bdata_def_cred

label var ymod_nd   "Model non-default (output)"
label var ymod_def  "Model default-linked (output)"
label var lmod_nd   "Model non-default (credit)"
label var lmod_def  "Model default-linked (credit)"
label var bdata_nd  "Empirical non-default (output)"
label var bdata_def "Empirical default-linked (output)"
save "$clean/cal_transmission.dta", replace

* ── Figure 1: Output — model vs. data ─────────────────────────────────────
twoway ///
    (connected bdata_nd  horizon, lcolor("0 84 166")  mcolor("0 84 166")  msymbol(circle)) ///
    (connected ymod_nd   horizon, lcolor("0 84 166")  lpattern(dash)      msymbol(none)) ///
    (connected bdata_def horizon, lcolor("157 36 73") mcolor("157 36 73") msymbol(square)) ///
    (connected ymod_def  horizon, lcolor("157 36 73") lpattern(dash)      msymbol(none)), ///
    yline(0, lpattern(dash) lcolor(gs10)) ///
    xlabel(0(1)4) ylabel(, format(%4.1f)) ///
    xtitle("Years after onset") ytitle("Cumulative output response (pp)") ///
    title("Model vs. Data: Output IRFs", size(medium)) ///
    subtitle("Solid = LP estimate; dashed = calibrated model", size(small)) ///
    legend(order(1 "Data: non-default" 2 "Model: non-default" ///
                 3 "Data: default" 4 "Model: default") size(vsmall) rows(2)) ///
    note("ND params {xi,phi,Phi_N} SMM-fit; default path out-of-sample.", size(vsmall)) ///
    graphregion(color(white))
graph export "$figs/fig_model_vs_data_output.pdf", replace
graph export "$figs/fig_model_vs_data_output.png", replace width(1200)

* ── Figure 2: Credit — model vs. data (if credit IRFs available) ───────────
capture confirm variable bdata_nd_cred
if _rc == 0 {
    twoway ///
        (connected bdata_nd_cred  horizon, lcolor("0 84 166")  mcolor("0 84 166")  msymbol(circle)) ///
        (connected lmod_nd        horizon, lcolor("0 84 166")  lpattern(dash)      msymbol(none)) ///
        (connected bdata_def_cred horizon, lcolor("157 36 73") mcolor("157 36 73") msymbol(square)) ///
        (connected lmod_def       horizon, lcolor("157 36 73") lpattern(dash)      msymbol(none)), ///
        yline(0, lpattern(dash) lcolor(gs10)) ///
        xlabel(0(1)4) ylabel(, format(%4.1f)) ///
        xtitle("Years after onset") ytitle("Cumulative private credit/GDP response (pp)") ///
        title("Model vs. Data: Private Credit IRFs", size(medium)) ///
        subtitle("Solid = LP estimate; dashed = model (ND: leverage amplifier; DEF: ≈0)", size(small)) ///
        legend(order(1 "Data: non-default" 2 "Model: non-default" ///
                     3 "Data: default" 4 "Model: default") size(vsmall) rows(2)) ///
        note("ND: credit via balance-sheet (IR.4); DEF: credit/GDP≈0 (GDP-denominator effect).", size(vsmall)) ///
        graphregion(color(white))
    graph export "$figs/fig_model_vs_data_credit.pdf", replace
    graph export "$figs/fig_model_vs_data_credit.png", replace width(1200)
    di as result "Credit figure saved: fig_model_vs_data_credit.pdf"
}
else {
    di as text "Credit IRF files not found — credit figure skipped."
    di as text "Run 03_lp_resolution.do with credit outcome to generate irf_nd_credit.dta / irf_def_credit.dta"
}

di as result _n "════════════════════════════════════════════════════════"
di as result "TRANSMISSION BLOCK COMPLETE."
di as result "Fitted: xi=" %4.3f xi_fit "  phi=" %4.3f phi_fit "  Phi_N=" %4.3f PhiN_fit
di as result "Figure: fig_model_vs_data.pdf"
di as result "════════════════════════════════════════════════════════"
