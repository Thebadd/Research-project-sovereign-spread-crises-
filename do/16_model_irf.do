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
beta   = st_numscalar("beta")
alpha  = st_numscalar("alpha")
delta  = st_numscalar("delta")
s_ss   = st_numscalar("s_ss")
bBN    = st_numscalar("bB_gdp")/st_numscalar("N_gdp")   // bB_ss / N_ss
credit_gdp = st_numscalar("credit_gdp")
N_gdp  = st_numscalar("N_gdp")
lambda = st_numscalar("lambda")
Phi_O  = st_numscalar("Phi_O")
gamma  = st_numscalar("gamma")          // sovereign ceiling, set externally

// spread path (decimal), pick horizons 0..4 (rows where horizon>=0)
SP = st_matrix("SPATH")                 // col1 horizon, col2 s_model
shat = J(5,1,.)
for (h=0; h<=4; h++) {
    idx = selectindex(SP[.,1]:==h)
    shat[h+1] = SP[idx,2] - s_ss        // deviation from SS spread
}

// empirical non-default output IRF, h=0..4
ND = st_matrix("IRF_ND")
bnd = J(5,1,.)
for (h=0;h<=4;h++) { idx=selectindex(ND[.,1]:==h); bnd[h+1]=ND[idx,2] }

// ─── Transmission IRF given deep params ────────────────────────────────────
// returns 5x1 model output IRF (in pp) for h=0..4
function trans_irf(xi, phi, Phi_N, shat, gamma, alpha, delta, beta, bBN, Phi_O) {
    // steady-state objects depending on (xi, phi)
    RL_ss = 1/beta + phi*bBN
    nexp  = (1-alpha)/alpha
    eps_p = xi*nexp*RL_ss/(1 + xi*(RL_ss-1))
    eta   = 1/((1-alpha)*(RL_ss-(1-delta)))
    Omega = phi*bBN
    Bsens = bBN/(1+st_numscalar("s_ss"))^2

    H = 5
    nw  = J(H,1,0)
    kap = J(H,1,0)
    yhat= J(H,1,0)

    // h=0 (index 1)
    nw[1]  = -Bsens*shat[1]
    kap[1] = 0
    yhat[1]= alpha*kap[1] - eps_p*(gamma*shat[1] - Omega*nw[1])
    // h>=1
    for (h=2; h<=H; h++) {
        nw[h]  = Phi_N*nw[h-1] - Bsens*shat[h] - Phi_O*gamma*shat[h-1]
        kap[h] = (1-delta)*kap[h-1] - eta*(gamma*shat[h-1] - Omega*nw[h-1])
        yhat[h]= alpha*kap[h] - eps_p*(gamma*shat[h] - Omega*nw[h])
    }
    return(yhat:*100)               // -> percentage points
}

// ─── SMM grid search over {xi, phi, Phi_N} ────────────────────────────────
xig = range(0.20,0.90,0.05)
phg = range(0.01,0.15,0.01)
png = range(0.30,0.95,0.05)

best = 1e12; bx=.; bp=.; bn=.
for (i=1;i<=rows(xig);i++) {
 for (j=1;j<=rows(phg);j++) {
  for (k=1;k<=rows(png);k++) {
     yh = trans_irf(xig[i], phg[j], png[k], shat, gamma, alpha, delta, beta, bBN, Phi_O)
     sse = sum((yh - bnd):^2)
     if (sse < best) { best=sse; bx=xig[i]; bp=phg[j]; bn=png[k] }
  }
 }
}
printf("\n=== SMM FIT (non-default output IRF) ===\n")
printf("  xi*    = %5.3f\n", bx)
printf("  phi*   = %5.3f\n", bp)
printf("  Phi_N* = %5.3f\n", bn)
printf("  SSE    = %8.4f\n", best)

yfit = trans_irf(bx, bp, bn, shat, gamma, alpha, delta, beta, bBN, Phi_O)
printf("\n  h   data(nd)   model\n")
for (h=0;h<=4;h++) printf(" %2.0f   %7.3f   %7.3f\n", h, bnd[h+1], yfit[h+1])

st_matrix("YFIT", yfit)
st_numscalar("xi_fit", bx)
st_numscalar("phi_fit", bp)
st_numscalar("PhiN_fit", bn)

// ═══════════════════════════════════════════════════════════════════════════
//  DEFAULT PATH (validation, out-of-sample)
//  IR.5 corrected: under autarky I=0 so capital depletes by depreciation:
//     k_h^def = h*ln(1-delta)        (pure depreciation, no investment)
//  output: y_h^def = alpha*k_h^def - eps_p*dRL_aut
//  dRL_aut (autarky lending-rate wedge) is the ONE default-specific scalar;
//  it is pinned to the SINGLE h=0 default moment, then h=1..4 are checked
//  WITHOUT further tuning -> genuine out-of-sample validation.
// ═══════════════════════════════════════════════════════════════════════════
DEF = st_matrix("IRF_DEF")
bdef = J(5,1,.)
for (h=0;h<=4;h++) { idx=selectindex(DEF[.,1]:==h); bdef[h+1]=DEF[idx,2] }

// recompute eps_p at fitted params
RL_ss = 1/beta + bp*bBN
nexp  = (1-alpha)/alpha
eps_p = bx*nexp*RL_ss/(1+bx*(RL_ss-1))

// pin dRL_aut to h=0 default: bdef[1] = alpha*0 - eps_p*dRL_aut*100
dRL_aut = -bdef[1]/(eps_p*100)
printf("\n=== DEFAULT PATH VALIDATION ===\n")
printf("  Autarky lending-rate wedge dRL_aut (pinned to h=0) = %6.4f\n", dRL_aut)

ydef = J(5,1,.)
printf("\n  h   data(def)  model(def)   [h>=1 = out-of-sample]\n")
for (h=0;h<=4;h++) {
    kdef = h*ln(1-delta)
    ydef[h+1] = (alpha*kdef - eps_p*dRL_aut)*100
    tag = (h==0 ? "(pinned)" : "(oos)")
    printf(" %2.0f   %7.3f   %7.3f    %s\n", h, bdef[h+1], ydef[h+1], tag)
}
st_matrix("YDEF", ydef)

end

* ── Assemble model-vs-data dataset ────────────────────────────────────────
clear
set obs 5
gen horizon = _n - 1
svmat YFIT, names(ymod_nd)
svmat YDEF, names(ymod_def)
rename ymod_nd1  ymod_nd
rename ymod_def1 ymod_def

* merge empirical
merge 1:1 horizon using "$clean/irf_nd.dta", keepusing(b) nogen keep(master match)
rename b bdata_nd
merge 1:1 horizon using "$clean/irf_def.dta", keepusing(b) nogen keep(master match)
rename b bdata_def

label var ymod_nd   "Model non-default"
label var ymod_def  "Model default-linked"
label var bdata_nd  "Empirical non-default"
label var bdata_def "Empirical default-linked"
save "$clean/cal_transmission.dta", replace

* ── Overlay figure: model vs. data ─────────────────────────────────────────
twoway ///
    (connected bdata_nd  horizon, lcolor("0 84 166")   mcolor("0 84 166")   msymbol(circle)) ///
    (connected ymod_nd   horizon, lcolor("0 84 166")   lpattern(dash) msymbol(none)) ///
    (connected bdata_def horizon, lcolor("157 36 73")  mcolor("157 36 73")  msymbol(square)) ///
    (connected ymod_def  horizon, lcolor("157 36 73")  lpattern(dash) msymbol(none)), ///
    yline(0, lpattern(dash) lcolor(gs10)) ///
    xlabel(0(1)4) ylabel(, format(%4.1f)) ///
    xtitle("Years after onset") ytitle("Cumulative output response (pp)") ///
    title("Model vs. Empirical IRFs", size(medium)) ///
    subtitle("Solid = data (LP); dashed = calibrated model", size(small)) ///
    legend(order(1 "Data: non-default" 2 "Model: non-default" ///
                 3 "Data: default" 4 "Model: default") size(vsmall) rows(2)) ///
    note("Non-default params {xi,phi,Phi_N} SMM-fit to data; default path validated out-of-sample.", size(vsmall)) ///
    graphregion(color(white))
graph export "$figs/fig_model_vs_data.pdf", replace
graph export "$figs/fig_model_vs_data.png", replace width(1200)

di as result _n "════════════════════════════════════════════════════════"
di as result "TRANSMISSION BLOCK COMPLETE."
di as result "Fitted: xi=" %4.3f xi_fit "  phi=" %4.3f phi_fit "  Phi_N=" %4.3f PhiN_fit
di as result "Figure: fig_model_vs_data.pdf"
di as result "════════════════════════════════════════════════════════"
