#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import json

import numpy as np
import pandas as pd

ROOT = Path(".").resolve()

OUT = ROOT / (
    "runs/phase2/day051_f46ae_engineered_3x500ps_final"
)

AD = ROOT / (
    "runs/phase2/day050_f46ad_engineered_500ps_extensions"
)

F46R = ROOT / (
    "runs/phase2/day048_f46q_collective_200ps_extension/"
    "F46R_collective_200ps/analyze_F46R.py"
)

BULK_R1_NPZ = ROOT / (
    "runs/phase2/day049_f46t_bulk_500ps_convergence/"
    "F46T_bulk_collective_series_0_500.npz"
)

BULK_REP_ROOT = ROOT / (
    "runs/phase2/day049_f46w_bulk_independent_replicas"
)

# ------------------------------------------------------------------
# Validated F46R method
# ------------------------------------------------------------------

spec = importlib.util.spec_from_file_location(
    "f46r",
    F46R
)

f46r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f46r)

print("="*88)
print("F46AE — FINAL ENGINEERED 3×500 ps CONVERGENCE")
print("="*88)
print("F46R_VALIDATED_METHOD_IMPORT=PASS")


# ------------------------------------------------------------------
# Load helpers
# ------------------------------------------------------------------

def frames_to_arrays(t, frames, nw):

    t = np.asarray(t, dtype=float)

    M = np.asarray(
        [f["M_e_nm"] for f in frames],
        dtype=float
    )

    boxes = np.asarray(
        [f["box_nm"] for f in frames],
        dtype=float
    )

    return t, M, boxes, int(nw)


def load_xtc(tpr, xtc):

    t,frames,nw = f46r.load_segment(
        tpr,
        xtc
    )

    return frames_to_arrays(
        t,frames,nw
    )


engineered = {}

for rep in ["R1","R2","R3"]:

    engineered[rep] = load_xtc(
        AD/rep/f"{rep}_engineered_500ps.tpr",
        OUT/rep/f"{rep}_engineered_0_500ps.xtc",
    )


# ------------------------------------------------------------------
# Strict input validation
# ------------------------------------------------------------------

for rep,(t,M,boxes,nw) in engineered.items():

    dt=np.diff(t)

    gate=(
        nw == 16634
        and len(t) == 1001
        and abs(t[0]) < 1e-8
        and abs(t[-1]-500.0) < 1e-8
        and len(np.unique(np.round(t,6))) == 1001
        and np.allclose(
            dt,
            0.5,
            atol=1e-6,
            rtol=0
        )
    )

    print(
        f"ENGINEERED_{rep}_500PS_INPUT_GATE="
        + ("PASS" if gate else "FAIL")
    )

    if not gate:
        raise RuntimeError(
            f"{rep}: invalid 0-500 ps input"
        )

print(
    "F46AE_ENGINEERED_INPUT_GATE=PASS"
)


# ------------------------------------------------------------------
# Windows
#
# Cumulative windows assess convergence.
# Late windows distinguish persistent behavior from early transient.
# ------------------------------------------------------------------

ENG_WINDOWS = [
    ("0_300",0.0,300.0),
    ("0_400",0.0,400.0),
    ("0_500",0.0,500.0),
    ("100_500",100.0,500.0),
    ("200_500",200.0,500.0),
    ("300_500",300.0,500.0),
    ("400_500",400.0,500.0),
]

MATCHED_BULK_WINDOWS = [
    ("0_300",0.0,300.0),
    ("100_300",100.0,300.0),
    ("150_300",150.0,300.0),
    ("200_300",200.0,300.0),
]


# ------------------------------------------------------------------
# Tensor helpers
# ------------------------------------------------------------------

def tensor_result(M,boxes):

    T,V = f46r.fluctuation_tensor(
        M,boxes
    )

    met = f46r.tensor_metrics(T)

    eigvals,eigvecs = np.linalg.eigh(T)

    order=np.argsort(eigvals)[::-1]

    eigvals=eigvals[order]
    eigvecs=eigvecs[:,order]

    return T,V,met,eigvals,eigvecs


tensor_rows=[]
tensor_cache={}

for rep,(t,M,boxes,nw) in engineered.items():

    tensor_cache[rep]={}

    for name,start,end in ENG_WINDOWS:

        mask=(
            (t>=start-1e-8)
            & (t<=end+1e-8)
        )

        T,V,met,eigvals,eigvecs = tensor_result(
            M[mask],
            boxes[mask]
        )

        tensor_cache[rep][name]={
            "T":T,
            "eigvals":eigvals,
            "eigvecs":eigvecs,
        }

        tensor_rows.append({
            "replica":rep,
            "window":name,
            "start_ps":start,
            "end_ps":end,
            "n_frames":int(mask.sum()),
            "volume_nm3":V,
            **met,
            "eig1":float(eigvals[0]),
            "eig2":float(eigvals[1]),
            "eig3":float(eigvals[2]),
        })


tensor_df=pd.DataFrame(
    tensor_rows
)


# ------------------------------------------------------------------
# Equal-weight ensemble tensors
# ------------------------------------------------------------------

ensemble_rows=[]
ensemble_cache={}

for name,start,end in ENG_WINDOWS:

    Ts=[
        tensor_cache[r][name]["T"]
        for r in ["R1","R2","R3"]
    ]

    Tmean=sum(Ts)/3.0

    met=f46r.tensor_metrics(
        Tmean
    )

    eigvals,eigvecs=np.linalg.eigh(
        Tmean
    )

    order=np.argsort(eigvals)[::-1]

    eigvals=eigvals[order]
    eigvecs=eigvecs[:,order]

    ensemble_cache[name]={
        "T":Tmean,
        "eigvals":eigvals,
        "eigvecs":eigvecs,
    }

    ensemble_rows.append({
        "window":name,
        "start_ps":start,
        "end_ps":end,
        **met,
        "eig1":float(eigvals[0]),
        "eig2":float(eigvals[1]),
        "eig3":float(eigvals[2]),
        "gap12_fraction":
            float(
                abs(eigvals[0]-eigvals[1])
                / max(abs(eigvals[0]),1e-15)
            ),
        "gap23_fraction":
            float(
                abs(eigvals[1]-eigvals[2])
                / max(abs(eigvals[1]),1e-15)
            ),
    })


ensemble_df=pd.DataFrame(
    ensemble_rows
)


# ------------------------------------------------------------------
# Between-replica uncertainty
# ------------------------------------------------------------------

OBS = [
    "XX","YY","ZZ",
    "XY","XZ","YZ",
    "trace",
    "eig1","eig2","eig3",
]

stat_rows=[]

for window,_,_ in ENG_WINDOWS:

    sub=tensor_df[
        tensor_df.window==window
    ]

    for obs in OBS:

        x=sub[obs].to_numpy(
            dtype=float
        )

        mean=float(np.mean(x))
        sd=float(np.std(x,ddof=1))
        sem=float(sd/np.sqrt(3))

        stat_rows.append({
            "window":window,
            "observable":obs,
            "mean":mean,
            "sd":sd,
            "sem":sem,
            "cv":
                (
                    sd/abs(mean)
                    if abs(mean)>1e-15
                    else np.nan
                ),
        })


stats_df=pd.DataFrame(
    stat_rows
)


# ------------------------------------------------------------------
# PRIMARY convergence criterion
#
# No isotropy gate for engineered.
#
# Compare 0-400 -> 0-500 change against 500-ps
# between-replica SD.
#
# Frozen uncertainty-based principle:
# temporal change must not exceed replica uncertainty.
# ------------------------------------------------------------------

PRIMARY=[
    "trace",
    "eig1",
    "eig2",
    "eig3",
]

conv_rows=[]

for obs in PRIMARY:

    v400=float(
        ensemble_df[
            ensemble_df.window=="0_400"
        ][obs].iloc[0]
    )

    v500=float(
        ensemble_df[
            ensemble_df.window=="0_500"
        ][obs].iloc[0]
    )

    sd500=float(
        stats_df[
            (stats_df.window=="0_500")
            & (stats_df.observable==obs)
        ].sd.iloc[0]
    )

    delta=abs(v500-v400)

    gate=(
        delta <= sd500
        if sd500>0
        else False
    )

    conv_rows.append({
        "observable":obs,
        "ensemble_0_400":v400,
        "ensemble_0_500":v500,
        "absolute_change":delta,
        "relative_change_percent":
            (
                100*delta/abs(v400)
                if abs(v400)>1e-15
                else np.nan
            ),
        "between_replica_sd_0_500":sd500,
        "change_over_sd":
            (
                delta/sd500
                if sd500>0
                else np.nan
            ),
        "gate":gate,
    })


conv_df=pd.DataFrame(
    conv_rows
)

tensor_convergence_gate=bool(
    conv_df.gate.all()
)


# ------------------------------------------------------------------
# Effective sample size
# Same practical diagnostic as F46AC.
# ------------------------------------------------------------------

def acf_fft(x):

    x=np.asarray(x,dtype=float)
    x=x-np.mean(x)

    n=len(x)

    var=np.dot(x,x)/n

    if var<=0:
        return np.ones(n)

    nfft=1 << (2*n-1).bit_length()

    f=np.fft.rfft(
        x,
        n=nfft
    )

    ac=np.fft.irfft(
        f*np.conjugate(f),
        n=nfft
    )[:n]

    ac=ac/np.arange(
        n,0,-1,
        dtype=float
    )

    return ac/ac[0]


def tau_neff(x,dt):

    rho=acf_fft(x)

    idx=np.where(
        rho[1:]<=0
    )[0]

    if len(idx):
        stop=int(idx[0])+1
    else:
        stop=len(rho)-1

    tau=dt*(
        0.5
        + np.sum(
            rho[1:stop+1]
        )
    )

    tau=max(
        float(tau),
        dt/2
    )

    neff=len(x)*dt/(2*tau)

    return tau,float(neff)


neff_rows=[]

for rep,(t,M,boxes,nw) in engineered.items():

    for name,start,end in [
        ("0_500",0,500),
        ("200_500",200,500),
        ("300_500",300,500),
    ]:

        mask=(
            (t>=start-1e-8)
            & (t<=end+1e-8)
        )

        X=M[mask]

        dM=X-np.mean(
            X,
            axis=0,
            keepdims=True
        )

        q={
            "Mx2":dM[:,0]**2,
            "My2":dM[:,1]**2,
            "Mz2":dM[:,2]**2,
            "M2_trace":
                np.sum(
                    dM*dM,
                    axis=1
                ),
        }

        dt=float(
            np.median(
                np.diff(t[mask])
            )
        )

        for obs,x in q.items():

            tau,neff=tau_neff(
                x,dt
            )

            neff_rows.append({
                "replica":rep,
                "window":name,
                "observable":obs,
                "tau_int_ps":tau,
                "Neff":neff,
            })


neff_df=pd.DataFrame(
    neff_rows
)

min_neff_500=float(
    neff_df[
        neff_df.window=="0_500"
    ].Neff.min()
)

# Practical sampling diagnostic only, not theorem.
neff_gate = (
    min_neff_500 >= 30.0
)


# ------------------------------------------------------------------
# Collective correlations
# ------------------------------------------------------------------

corr_rows=[]

for rep,(t,M,boxes,nw) in engineered.items():

    for name,start,end in [
        ("0_400",0,400),
        ("0_500",0,500),
        ("200_500",200,500),
        ("300_500",300,500),
        ("400_500",400,500),
    ]:

        mask=(
            (t>=start-1e-8)
            & (t<=end+1e-8)
        )

        tt=t[mask]-t[mask][0]
        X=M[mask]

        C,origins=f46r.connected_corr(X)

        iz,zt=f46r.first_zero_integral(
            tt,C
        )

        corr_rows.append({
            "replica":rep,
            "window":name,
            "integral_5ps":
                f46r.integral_to(
                    tt,C,5
                ),
            "integral_10ps":
                f46r.integral_to(
                    tt,C,10
                ),
            "first_zero_integral_ps":
                iz,
            "first_zero_time_ps":
                zt,
        })


corr_df=pd.DataFrame(
    corr_rows
)


# ------------------------------------------------------------------
# Correlation convergence 0-400 -> 0-500
# ------------------------------------------------------------------

corr_conv_rows=[]

for obs in [
    "integral_5ps",
    "integral_10ps",
    "first_zero_integral_ps",
]:

    a=corr_df[
        corr_df.window=="0_400"
    ][obs].to_numpy(dtype=float)

    b=corr_df[
        corr_df.window=="0_500"
    ][obs].to_numpy(dtype=float)

    mean400=float(np.mean(a))
    mean500=float(np.mean(b))

    sd500=float(np.std(b,ddof=1))

    delta=abs(mean500-mean400)

    gate=(
        delta <= sd500
        if sd500>0
        else False
    )

    corr_conv_rows.append({
        "observable":obs,
        "mean_0_400":mean400,
        "mean_0_500":mean500,
        "absolute_change":delta,
        "sd_between_replica_0_500":sd500,
        "change_over_sd":
            (
                delta/sd500
                if sd500>0
                else np.nan
            ),
        "gate":gate,
    })


corr_conv_df=pd.DataFrame(
    corr_conv_rows
)

correlation_gate=bool(
    corr_conv_df.gate.all()
)


# ------------------------------------------------------------------
# Principal-axis reproducibility
# Direction is diagnostic only if eigenvalue is non-degenerate.
# ------------------------------------------------------------------

axis_rows=[]

for window in [
    "0_500",
    "200_500",
    "300_500",
]:

    Vmean=ensemble_cache[
        window
    ]["eigvecs"]

    evals=ensemble_cache[
        window
    ]["eigvals"]

    for k in range(3):

        if k==0:
            gap=abs(evals[0]-evals[1])/abs(evals[0])

        elif k==2:
            gap=abs(evals[1]-evals[2])/abs(evals[1])

        else:
            gap=min(
                abs(evals[0]-evals[1])/abs(evals[1]),
                abs(evals[1]-evals[2])/abs(evals[1]),
            )

        for rep in ["R1","R2","R3"]:

            V=tensor_cache[
                rep
            ][window]["eigvecs"]

            dot=float(
                np.clip(
                    abs(
                        np.dot(
                            V[:,k],
                            Vmean[:,k]
                        )
                    ),
                    0,1
                )
            )

            angle=float(
                np.degrees(
                    np.arccos(dot)
                )
            )

            axis_rows.append({
                "window":window,
                "eigen_index":k+1,
                "replica":rep,
                "ensemble_gap_fraction":
                    float(gap),
                "angle_deg":angle,
                "interpretation":
                    (
                        "INTERPRETABLE"
                        if gap>=0.10
                        else
                        "NEAR_DEGENERATE_DIAGNOSTIC_ONLY"
                    ),
            })


axis_df=pd.DataFrame(
    axis_rows
)


# ------------------------------------------------------------------
# Low-frequency collective spectral fractions, 0-1 THz.
# NOT excitonic spectral density.
# ------------------------------------------------------------------

BANDS=[
    ("0_0p05",0.00,0.05),
    ("0p05_0p20",0.05,0.20),
    ("0p20_0p50",0.20,0.50),
    ("0p50_1p00",0.50,1.00),
]

spectral_rows=[]

for rep,(t,M,boxes,nw) in engineered.items():

    for window,start,end in [
        ("0_500",0,500),
        ("200_500",200,500),
        ("300_500",300,500),
    ]:

        mask=(
            (t>=start-1e-8)
            & (t<=end+1e-8)
        )

        tt=t[mask]
        X=M[mask]

        dt=float(
            np.median(np.diff(tt))
        )

        X=X-np.mean(
            X,
            axis=0,
            keepdims=True
        )

        win=np.hanning(len(X))

        F=np.fft.rfft(
            X*win[:,None],
            axis=0
        )

        freq=np.fft.rfftfreq(
            len(X),
            d=dt
        )

        P=np.sum(
            np.abs(F)**2,
            axis=1
        )

        total=float(
            np.trapezoid(
                P,freq
            )
        )

        for band,lo,hi in BANDS:

            m=(
                (freq>=lo)
                & (
                    freq<=hi
                    if hi==1.0
                    else freq<hi
                )
            )

            power=float(
                np.trapezoid(
                    P[m],
                    freq[m]
                )
            )

            spectral_rows.append({
                "replica":rep,
                "window":window,
                "band":band,
                "power_arbitrary":power,
                "fraction_0_1THz":
                    (
                        power/total
                        if total>0
                        else np.nan
                    ),
                "nyquist_THz":
                    float(freq[-1]),
            })


spectral_df=pd.DataFrame(
    spectral_rows
)


# ==================================================================
# MATCHED ENGINEERED-vs-BULK comparison
#
# ONLY <=300 ps because bulk R2/R3 are 300 ps.
# ==================================================================

bulk={}

d=np.load(
    BULK_R1_NPZ
)

bulk["R1"]=(
    np.asarray(d["time_ps"],dtype=float),
    np.asarray(d["M"],dtype=float),
    np.asarray(d["boxes"],dtype=float),
    16634,
)

for rep in ["R2","R3"]:

    bulk[rep]=load_xtc(
        BULK_REP_ROOT/rep/"production_300ps.tpr",
        BULK_REP_ROOT/rep/"production_300ps.xtc",
    )


comparison_rows=[]

for window,start,end in MATCHED_BULK_WINDOWS:

    eng_T=[]
    bulk_T=[]

    eng_rep_metrics=[]
    bulk_rep_metrics=[]

    for rep in ["R1","R2","R3"]:

        te,Me,be,nwe=engineered[rep]

        me=(
            (te>=start-1e-8)
            & (te<=end+1e-8)
        )

        Te,Ve,mete,eige,vece = tensor_result(
            Me[me],
            be[me]
        )

        eng_T.append(Te)

        eng_rep_metrics.append({
            **mete,
            "eig1":float(eige[0]),
            "eig2":float(eige[1]),
            "eig3":float(eige[2]),
        })


        tb,Mb,bb,nwb=bulk[rep]

        mb=(
            (tb>=start-1e-8)
            & (tb<=end+1e-8)
        )

        Tb,Vb,metb,eigb,vecb = tensor_result(
            Mb[mb],
            bb[mb]
        )

        bulk_T.append(Tb)

        bulk_rep_metrics.append({
            **metb,
            "eig1":float(eigb[0]),
            "eig2":float(eigb[1]),
            "eig3":float(eigb[2]),
        })


    E=sum(eng_T)/3.0
    B=sum(bulk_T)/3.0

    Emet=f46r.tensor_metrics(E)
    Bmet=f46r.tensor_metrics(B)

    Eeig=np.sort(
        np.linalg.eigvalsh(E)
    )[::-1]

    Beig=np.sort(
        np.linalg.eigvalsh(B)
    )[::-1]

    Eensemble={
        **Emet,
        "eig1":float(Eeig[0]),
        "eig2":float(Eeig[1]),
        "eig3":float(Eeig[2]),
    }

    Bensemble={
        **Bmet,
        "eig1":float(Beig[0]),
        "eig2":float(Beig[1]),
        "eig3":float(Beig[2]),
    }

    for obs in [
        "XX","YY","ZZ",
        "trace",
        "eig1","eig2","eig3",
    ]:

        e=float(Eensemble[obs])
        b=float(Bensemble[obs])

        er=np.asarray(
            [x[obs] for x in eng_rep_metrics],
            dtype=float
        )

        br=np.asarray(
            [x[obs] for x in bulk_rep_metrics],
            dtype=float
        )

        delta_rep=er-br

        same_sign=max(
            int(np.sum(delta_rep>0)),
            int(np.sum(delta_rep<0))
        )

        sd_comb=np.sqrt(
            np.std(er,ddof=1)**2
            + np.std(br,ddof=1)**2
        )

        comparison_rows.append({
            "window":window,
            "observable":obs,
            "engineered":e,
            "bulk":b,
            "difference":e-b,
            "percent_difference":
                (
                    100*(e-b)/abs(b)
                    if abs(b)>1e-15
                    else np.nan
                ),
            "same_sign_replica_count":
                same_sign,
            "combined_between_replica_sd":
                float(sd_comb),
            "difference_over_combined_sd":
                (
                    abs(e-b)/sd_comb
                    if sd_comb>0
                    else np.nan
                ),
        })


comparison_df=pd.DataFrame(
    comparison_rows
)


# ------------------------------------------------------------------
# Persistent / transient classification
#
# Requires late matched 200-300 behavior for persistent claim.
# ------------------------------------------------------------------

class_rows=[]

for obs in [
    "XX","YY","ZZ",
    "trace",
    "eig1","eig2","eig3",
]:

    full=comparison_df[
        (comparison_df.window=="0_300")
        & (comparison_df.observable==obs)
    ].iloc[0]

    late=comparison_df[
        (comparison_df.window=="200_300")
        & (comparison_df.observable==obs)
    ].iloc[0]

    same=int(
        late.same_sign_replica_count
    )

    sep=float(
        late.difference_over_combined_sd
    )

    full_pct=float(
        full.percent_difference
    )

    late_pct=float(
        late.percent_difference
    )

    if (
        same==3
        and np.isfinite(sep)
        and sep>=1.0
    ):
        classification=(
            "PERSISTENT_REPLICA_REPRODUCIBLE"
        )

    elif (
        abs(full_pct)>=10.0
        and abs(late_pct)
            < 0.5*abs(full_pct)
    ):
        classification=(
            "PREDOMINANTLY_TRANSIENT"
        )

    else:
        classification=(
            "UNRESOLVED_OR_WEAK"
        )

    class_rows.append({
        "observable":obs,
        "full_0_300_percent_difference":
            full_pct,
        "late_200_300_percent_difference":
            late_pct,
        "late_same_sign_replica_count":
            same,
        "late_difference_over_combined_sd":
            sep,
        "classification":
            classification,
    })


class_df=pd.DataFrame(
    class_rows
)


# ------------------------------------------------------------------
# Final decision
# ------------------------------------------------------------------

engineered_convergence = (
    tensor_convergence_gate
    and neff_gate
)

if engineered_convergence:

    decision=(
        "ENGINEERED_COLLECTIVE_RESPONSE_CONVERGED"
    )

else:

    decision=(
        "ENGINEERED_COLLECTIVE_RESPONSE_"
        "CONVERGENCE_NOT_ESTABLISHED"
    )


# ------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------

tensor_df.to_csv(
    OUT/"F46AE_replica_tensor_windows.csv",
    index=False
)

ensemble_df.to_csv(
    OUT/"F46AE_ensemble_tensor_windows.csv",
    index=False
)

stats_df.to_csv(
    OUT/"F46AE_between_replica_statistics.csv",
    index=False
)

conv_df.to_csv(
    OUT/"F46AE_primary_tensor_convergence.csv",
    index=False
)

neff_df.to_csv(
    OUT/"F46AE_effective_sample_size.csv",
    index=False
)

corr_df.to_csv(
    OUT/"F46AE_collective_correlations.csv",
    index=False
)

corr_conv_df.to_csv(
    OUT/"F46AE_correlation_convergence.csv",
    index=False
)

axis_df.to_csv(
    OUT/"F46AE_principal_axis_reproducibility.csv",
    index=False
)

spectral_df.to_csv(
    OUT/"F46AE_low_frequency_collective_spectrum.csv",
    index=False
)

comparison_df.to_csv(
    OUT/"F46AE_matched_engineered_vs_bulk.csv",
    index=False
)

class_df.to_csv(
    OUT/"F46AE_persistent_transient_classification.csv",
    index=False
)


summary={
    "engineered_dataset":{
        "replicas":3,
        "duration_ps":500,
        "dt_ps":0.5,
        "n_frames_each":1001,
    },
    "primary_tensor_convergence_gate":
        tensor_convergence_gate,
    "min_full500_Neff":
        min_neff_500,
    "Neff_practical_gate":
        neff_gate,
    "correlation_support_gate":
        correlation_gate,
    "engineered_convergence":
        engineered_convergence,
    "decision":
        decision,
    "comparison_boundary":
        (
            "Engineered-vs-bulk comparison is restricted "
            "to exactly matched <=300 ps windows because "
            "the converged three-replica bulk reference "
            "contains 300 ps per replica."
        ),
    "interpretation_boundaries":[
        (
            "No engineered-system isotropy criterion "
            "is imposed."
        ),
        (
            "Collective-polarization response is an "
            "environmental mechanism-level observable."
        ),
        (
            "No quantum coherence lifetime is inferred."
        ),
        (
            "0.5 ps sampling limits collective spectral "
            "interpretation to <=1 THz."
        ),
        (
            "The reported spectrum is not an excitonic "
            "spectral density."
        ),
        (
            "TIP4P/2005 lacks explicit electronic "
            "polarization."
        ),
        (
            "The frozen scaffold contains no dynamical "
            "scaffold polarization."
        ),
    ],
}


(
    OUT/"F46AE_FINAL_SUMMARY.json"
).write_text(
    json.dumps(
        summary,
        indent=2
    )+"\n"
)


# ==================================================================
# Terminal report
# ==================================================================

print()
print("="*88)
print("ENGINEERED INDIVIDUAL TENSORS — 0-500 ps")
print("="*88)

print(
    tensor_df[
        tensor_df.window=="0_500"
    ][
        [
            "replica",
            "XX","YY","ZZ",
            "trace",
            "eig1","eig2","eig3",
        ]
    ].to_string(index=False)
)


print()
print("="*88)
print("ENGINEERED EQUAL-WEIGHT ENSEMBLE — CONVERGENCE WINDOWS")
print("="*88)

print(
    ensemble_df[
        ensemble_df.window.isin(
            [
                "0_300",
                "0_400",
                "0_500",
                "200_500",
                "300_500",
                "400_500",
            ]
        )
    ][
        [
            "window",
            "XX","YY","ZZ",
            "trace",
            "eig1","eig2","eig3",
            "diag_cv",
        ]
    ].to_string(index=False)
)


print()
print("="*88)
print("PRIMARY TENSOR CONVERGENCE — 0-400 -> 0-500")
print("="*88)

print(
    conv_df.to_string(
        index=False
    )
)

print(
    "F46AE_PRIMARY_TENSOR_CONVERGENCE_GATE="
    + (
        "PASS"
        if tensor_convergence_gate
        else "FAIL"
    )
)


print()
print("="*88)
print("FULL 0-500 ps EFFECTIVE SAMPLE SIZE")
print("="*88)

print(
    neff_df[
        neff_df.window=="0_500"
    ].to_string(index=False)
)

print(
    f"MIN_FULL500_NEFF={min_neff_500}"
)

print(
    "F46AE_NEFF_GATE="
    + (
        "PASS"
        if neff_gate
        else "FAIL"
    )
)


print()
print("="*88)
print("COLLECTIVE-CORRELATION CONVERGENCE")
print("="*88)

print(
    corr_conv_df.to_string(
        index=False
    )
)

print(
    "F46AE_CORRELATION_SUPPORT_GATE="
    + (
        "PASS"
        if correlation_gate
        else "FAIL"
    )
)


print()
print("="*88)
print("PRINCIPAL AXES — ENGINEERED 0-500 ps")
print("="*88)

print(
    axis_df[
        axis_df.window=="0_500"
    ].to_string(index=False)
)


print()
print("="*88)
print("MATCHED ENGINEERED vs BULK — 0-300 ps")
print("="*88)

print(
    comparison_df[
        comparison_df.window=="0_300"
    ][
        [
            "observable",
            "engineered",
            "bulk",
            "percent_difference",
            "same_sign_replica_count",
            "difference_over_combined_sd",
        ]
    ].to_string(index=False)
)


print()
print("="*88)
print("MATCHED ENGINEERED vs BULK — LATE 200-300 ps")
print("="*88)

print(
    comparison_df[
        comparison_df.window=="200_300"
    ][
        [
            "observable",
            "engineered",
            "bulk",
            "percent_difference",
            "same_sign_replica_count",
            "difference_over_combined_sd",
        ]
    ].to_string(index=False)
)


print()
print("="*88)
print("PERSISTENT / TRANSIENT CLASSIFICATION")
print("="*88)

print(
    class_df.to_string(
        index=False
    )
)


print()
print("="*88)
print("FINAL F46AE ADJUDICATION")
print("="*88)

print(
    "F46AE_ENGINEERED_CONVERGENCE_GATE="
    + (
        "PASS"
        if engineered_convergence
        else "FAIL"
    )
)

print(
    f"F46AE_DECISION={decision}"
)

print(
    "NO_COHERENCE_LIFETIME_INFERRED=YES"
)

print(
    "F46AE_FINAL_ANALYSIS_COMPLETE=YES"
)
