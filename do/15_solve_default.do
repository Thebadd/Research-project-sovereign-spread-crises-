/*===========================================================================
  15_SOLVE_DEFAULT.DO
  Nonlinear default block (Theoretical Appendix Section 3, eqs 18-25).
  Canonical Arellano (2008) sovereign-default model with:
    - Tauchen-discretized AR(1) income process (calibrated rho_y, sigma_y)
    - Endogenous default decision (eq 24)
    - Risk-neutral foreign investors with recovery theta (eq 18)
    - Re-entry probability mu after exclusion (eq 23)
    - Bond-price fixed point q(B',y) (eqs 18-19, 25)

  Solved by value-function iteration in Mata.

  Outputs:
    - cal_default_moments.dta : model default freq, mean/sd spread, debt/GDP
    - cal_spread_event.dta    : event-study spread path h=-2..4 around onsets
                                (feeds the transmission block in 16_model_irf)

  REQUIRES: 14_calibration.do run first (cal_params.dta).

  NOTE: This is the genuinely nonlinear object — the default choice and the
        endogenous spread. The banking/output transmission is layered on top
        in 16_model_irf.do using the log-linear IRFs (IR.1-IR.5). This
        two-block design follows Bocola (2016): a sovereign-risk block plus
        a transmission block.
===========================================================================*/

use "$clean/cal_params.dta", clear

* Pull calibration into Mata-visible scalars
scalar beta   = beta[1]
scalar betap  = betap[1]
scalar sigma  = sigma[1]
scalar Rstar  = Rstar[1]
scalar theta  = theta[1]
scalar mu     = mu[1]
scalar rho_y  = rho_y[1]
scalar sigma_eps = sigma_eps[1]
scalar s_thr  = 0.10               // crisis threshold: 1000 bps = 0.10

* Default cost parameter (Arellano asymmetric cost): y_def = min(y, h_def*E[y]).
* Tunable to match default frequency (~2-4% annual in EM data).
scalar h_def  = 0.965

mata:
mata clear

// ─── Parameters from Stata scalars ───────────────────────────────────────
beta   = st_numscalar("beta")
betap  = st_numscalar("betap")
sigma  = st_numscalar("sigma")
Rstar  = st_numscalar("Rstar")
theta  = st_numscalar("theta")
mu     = st_numscalar("mu")
rho_y  = st_numscalar("rho_y")
sig_e  = st_numscalar("sigma_eps")
h_def  = st_numscalar("h_def")
s_thr  = st_numscalar("s_thr")

// ─── CRRA utility ─────────────────────────────────────────────────────────
function util(x, sigma) {
    // floor consumption to avoid log/neg-power blowups
    c = (x :< 1e-8) :* 1e-8 + (x :>= 1e-8) :* x
    if (sigma == 1) return(ln(c))
    return( (c:^(1-sigma) :- 1) :/ (1-sigma) )
}

// ─── Tauchen (1986) discretization of AR(1): log y' = rho*log y + eps ─────
function tauchen(rho, sig_e, N, m, ygrid, P) {
    sig_y = sig_e / sqrt(1 - rho^2)
    ymax  =  m*sig_y
    ymin  = -m*sig_y
    step  = (ymax - ymin)/(N-1)
    g = ymin :+ step:*(0::(N-1))           // grid in logs (Nx1)
    Pm = J(N, N, .)
    for (i=1; i<=N; i++) {
        for (j=1; j<=N; j++) {
            zlo = (g[j] - rho*g[i] - step/2)/sig_e
            zhi = (g[j] - rho*g[i] + step/2)/sig_e
            if (j==1)      Pm[i,j] = normal(zhi)
            else if (j==N) Pm[i,j] = 1 - normal(zlo)
            else           Pm[i,j] = normal(zhi) - normal(zlo)
        }
    }
    st_matrix(ygrid, exp(g))               // levels
    st_matrix(P, Pm)
}

// ─── Grids ────────────────────────────────────────────────────────────────
Ny = 21
Nb = 200
tauchen(rho_y, sig_e, Ny, 3, "ygr", "Pmat")
y = st_matrix("ygr")'                      // 1 x Ny (levels)
P = st_matrix("Pmat")                      // Ny x Ny

Ey = mean(y')                              // scalar E[y]
ydef = (y :< h_def*Ey) :* y + (y :>= h_def*Ey) :* (h_def*Ey)   // 1 x Ny autarky income

// Debt grid: B>0 = amount owed (face value), one-period discount bonds
Bmax = 1.5 * Ey
B = (0::(Nb-1)) :/ (Nb-1) :* Bmax          // Nb x 1

// ─── Initialize ───────────────────────────────────────────────────────────
VR = J(Nb, Ny, 0)        // value of repayment
VD = J(1,  Ny, 0)        // value of default (indep. of B)
V  = J(Nb, Ny, 0)        // V = max(VR, VD)
q  = J(Nb, Ny, 1/(1+Rstar))   // bond price q(B',y): rows = B' choice, cols = y today
defp = J(Nb, Ny, 0)      // default policy d(B,y)
bpol = J(Nb, Ny, 1)      // B' policy index

tol = 1e-6
maxit = 2000
disp = 1e9

for (it=1; it<=maxit; it++) {

    Vold = V
    qold = q

    // ----- Value of default (eq 23): autarky + re-entry w.p. mu to B=0 -----
    EVB0 = (V[1,.] * P')        // continuation if re-enter at B'=0: 1 x Ny  (V at B index 1)
    EVD  = (VD * P')            // continuation if still excluded:    1 x Ny
    VDnew = util(ydef, sigma) :+ betap :* ( mu:*EVB0 :+ (1-mu):*EVD )

    // ----- Value of repayment (eq 22): choose B' -----
    // resources given (B,y) and choice B': c = y - B + q(B',y)*B'
    EV = V * P'                 // Nb x Ny : E[V(B',y')|y] indexed by (B',y)
    VRnew = J(Nb, Ny, .)
    for (iy=1; iy<=Ny; iy++) {
        // consumption matrix: rows = current B, cols = candidate B'
        // c[ib, ibp] = y[iy] - B[ib] + q[ibp,iy]*B[ibp]
        qcol = q[.,iy]                          // Nb x 1
        issue = (qcol :* B)'                    // 1 x Nb : q(B')*B'
        cmat = (y[iy] :- B) :+ J(Nb,1,1)*issue  // Nb x Nb
        umat = util(cmat, sigma)                // Nb x Nb
        // add continuation E[V(B',y')|y] across candidate B' (cols)
        objmat = umat :+ betap :* J(Nb,1,1)*(EV[.,iy]')
        // maximize over B' (cols)
        maxv = rowmax(objmat)
        VRnew[.,iy] = maxv
        // record policy
        for (ib=1; ib<=Nb; ib++) {
            // argmax along the row
            rmax = max(objmat[ib,.])
            idx = selectindex(objmat[ib,.] :== rmax)[1]
            bpol[ib,iy] = idx
        }
    }

    // ----- Default decision (eq 24): default if VD > VR -----
    Vnew = J(Nb, Ny, .)
    for (iy=1; iy<=Ny; iy++) {
        for (ib=1; ib<=Nb; ib++) {
            if (VDnew[iy] > VRnew[ib,iy]) {
                Vnew[ib,iy] = VDnew[iy]
                defp[ib,iy] = 1
            }
            else {
                Vnew[ib,iy] = VRnew[ib,iy]
                defp[ib,iy] = 0
            }
        }
    }

    // ----- Bond price fixed point (eq 18): foreign investor zero profit -----
    // q(B',y) = (1/(1+r*)) * E_{y'|y}[ (1-d(B',y')) + d(B',y')*theta ]
    // payoff next period for holding B': 1 if repay, theta if default
    payoff = (1 :- defp) :+ defp:*theta        // Nb x Ny : payoff at (B', y')
    qnew = ( payoff * P' ) :/ (1+Rstar)        // Nb x Ny : expectation over y'|y

    // ----- Update & convergence -----
    V  = Vnew
    VR = VRnew
    VD = VDnew
    q  = 0.5:*qnew :+ 0.5:*qold                // dampened price update

    disp = max( (abs(V-Vold), abs(q-qold)) )
    if (mod(it,50)==0) printf("  iter %g   sup-norm = %g\n", it, disp)
    if (disp < tol) {
        printf("Converged in %g iterations (sup-norm = %g)\n", it, disp)
        break
    }
}

// ─── Store policy objects for simulation ────────────────────────────────────
st_matrix("M_q",    q)
st_matrix("M_def",  defp)
st_matrix("M_bpol", bpol)
st_matrix("M_y",    y')
st_matrix("M_P",    P)
st_matrix("M_B",    B)

// ═══════════════════════════════════════════════════════════════════════════
//  SIMULATION : long sample to compute model moments + crisis event study
// ═══════════════════════════════════════════════════════════════════════════
Tsim = 100000
burn = 1000
rseed(20260615)

// pre-compute cumulative transition for sampling y'
Pc = J(Ny, Ny, .)
for (i=1;i<=Ny;i++) Pc[i,.] = runningsum(P[i,.])

ysim = J(Tsim,1,.)
iy   = ceil(Ny/2)                  // start at median income
excl = 0                           // exclusion flag
ib   = 1                           // start at zero debt
ssim = J(Tsim,1,.)                 // spread
dsim = J(Tsim,1,.)                 // default event (new default)
bsim = J(Tsim,1,.)                 // debt/GDP
inmkt= J(Tsim,1,0)                 // in good standing?

for (t=1; t<=Tsim; t++) {
    // draw income
    u = runiform(1,1)
    iyn = selectindex(Pc[iy,.] :>= u)[1]
    ysim[t] = y[iyn]

    if (excl==0) {
        // in good standing: check default
        if (defp[ib,iyn]==1) {
            dsim[t]=1; excl=1; inmkt[t]=0; ssim[t]=.;
            ib=1;                 // debt wiped to grid zero on default
        }
        else {
            dsim[t]=0; inmkt[t]=1
            ibp = bpol[ib,iyn]
            // spread on newly issued bond: s = 1/q - 1 - r*
            qq = q[ibp,iyn]
            ssim[t] = (qq>1e-6 ? 1/qq - 1 - Rstar : .)
            bsim[t] = B[ibp]/ysim[t]
            ib = ibp
        }
    }
    else {
        // excluded: re-enter w.p. mu
        inmkt[t]=0; ssim[t]=.;
        if (runiform(1,1) < mu) { excl=0; ib=1; }
    }
    iy = iyn
}

// drop burn-in
ssim  = ssim[(burn+1)::Tsim]
dsim  = dsim[(burn+1)::Tsim]
bsim  = bsim[(burn+1)::Tsim]
inmkt = inmkt[(burn+1)::Tsim]
Ts = rows(ssim)

// ----- Model moments -----
def_freq = sum(dsim:==1)/Ts                            // annual default frequency
sm = select(ssim, inmkt:==1)                           // spreads when in market
mean_s = mean(sm)
sd_s   = sqrt(variance(sm))
bm = select(bsim, inmkt:==1 :& (bsim:!=.))
mean_b = mean(bm)

printf("\n=== MODEL MOMENTS ===\n")
printf("  Default frequency (annual) : %5.3f\n", def_freq)
printf("  Mean spread (decimal)      : %6.4f  (%4.0f bps)\n", mean_s, mean_s*10000)
printf("  SD spread (decimal)        : %6.4f  (%4.0f bps)\n", sd_s, sd_s*10000)
printf("  Mean debt/GDP (in market)  : %5.3f\n", mean_b)

st_numscalar("def_freq", def_freq)
st_numscalar("mean_s", mean_s)
st_numscalar("sd_s", sd_s)
st_numscalar("mean_b", mean_b)

// ----- Crisis event study: spread path around threshold crossings -----
// onset = first period spread crosses s_thr from below while in market
H = 4; Hpre = 2
acc = J(Hpre+H+1, 1, 0)      // sum of spreads by event-time
cnt = J(Hpre+H+1, 1, 0)      // count
for (t=(Hpre+1); t<=(Ts-H); t++) {
    if (inmkt[t]==1 & ssim[t]!=. & inmkt[t-1]==1 & ssim[t-1]!=.) {
        if (ssim[t] >= s_thr & ssim[t-1] < s_thr) {
            // event at t: collect window
            ok = 1
            for (k=-Hpre; k<=H; k++) if (ssim[t+k]==.) ok=0
            if (ok==1) {
                r = 0
                for (k=-Hpre; k<=H; k++) {
                    r = r+1
                    acc[r] = acc[r] + ssim[t+k]
                    cnt[r] = cnt[r] + 1
                }
            }
        }
    }
}
spath = acc :/ cnt           // mean spread by event time
hgrid = (-Hpre::H)
printf("\n=== CRISIS-EVENT SPREAD PATH (mean, decimal) ===\n")
for (r=1;r<=rows(spath);r++) printf("  h=%2.0f : %6.4f  (%4.0f bps)\n", hgrid[r], spath[r], spath[r]*10000)

st_matrix("M_spath", spath)
st_matrix("M_hgrid", hgrid)

end

* ── Export model moments ────────────────────────────────────────────────
clear
set obs 1
gen def_freq = def_freq
gen mean_s   = mean_s
gen sd_s     = sd_s
gen mean_b   = mean_b
save "$clean/cal_default_moments.dta", replace

* ── Export crisis-event spread path (feeds 16_model_irf.do) ──────────────
clear
svmat M_hgrid, names(horizon)
svmat M_spath, names(s_model)
rename horizon1 horizon
rename s_model1 s_model
label var s_model "Model-implied EMBIG spread (decimal)"
save "$clean/cal_spread_event.dta", replace

di as result _n "════════════════════════════════════════════════════════"
di as result "DEFAULT BLOCK SOLVED. Compare model vs. data moments:"
di as result "  Target (data): default freq ~ 0.03, mean spread ~ "  ///
    %4.0f (s_ss*10000) " bps"
di as result "  Tune h_def (default cost) and betap to match."
di as result "Spread event path saved -> cal_spread_event.dta (feeds 16)"
di as result "════════════════════════════════════════════════════════"
