"""Exploratory Cox survival model for CPTAC-PDA clinical data.

This model uses clinical predictors only. It does not reproduce the Jamalian
clinical-plus-radiomic model because radiomic feature values are not included
in the provided files. Outcomes:
  survival_days = Days to Death for Dead; Days to Last Follow Up for Alive
  event = 1 for Dead; 0 for Alive (right-censored)
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

def concordance_index(time, event, risk):
    time=np.asarray(time,float); event=np.asarray(event,int); risk=np.asarray(risk,float)
    concordant=discordant=tied=0
    for i in range(len(time)):
        if event[i] != 1:
            continue
        for j in range(len(time)):
            if time[i] < time[j]:
                if risk[i] > risk[j]: concordant += 1
                elif risk[i] < risk[j]: discordant += 1
                else: tied += 1
            elif i < j and time[i] == time[j] and event[j] == 1:
                tied += 1
    total=concordant+discordant+tied
    return float((concordant+0.5*tied)/total) if total else float("nan")

def fit_cox(X, time, event, penalizer=0.1):
    X=np.asarray(X,float); time=np.asarray(time,float); event=np.asarray(event,int)
    event_times=np.unique(time[event==1])
    n_events=max(int(event.sum()),1)
    def objective(beta):
        eta=np.clip(X @ beta, -30, 30)
        exp_eta=np.exp(eta)
        loss=0.0; grad=np.zeros(X.shape[1])
        for t in event_times:
            event_mask=(time==t)&(event==1)
            d=int(event_mask.sum())
            risk_mask=time>=t
            denom=max(exp_eta[risk_mask].sum(),1e-300)
            loss -= eta[event_mask].sum() - d*np.log(denom)
            grad -= X[event_mask].sum(axis=0) - d*(X[risk_mask]*exp_eta[risk_mask,None]).sum(axis=0)/denom
        loss=loss/n_events + 0.5*penalizer*np.dot(beta,beta)
        grad=grad/n_events + penalizer*beta
        return loss,grad
    res=minimize(objective,np.zeros(X.shape[1]),jac=True,method="BFGS")
    if not np.isfinite(res.fun):
        raise RuntimeError("Cox model optimization failed")
    return res.x

def make_features(df):
    x=pd.DataFrame(index=df.index)
    x["age_years"]=pd.to_numeric(df["age_years"],errors="coerce")
    x["female"]=(df["sex"].astype(str).str.lower()=="female").astype(float)
    x["high_grade"]=df["tumor_grade"].astype(str).str.upper().isin(["G3","G4"]).astype(float)
    return x

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data",default=str(ROOT/"data"/"cptac_pda_clinical_survival_clean.csv"))
    ap.add_argument("--folds",type=int,default=5)
    ap.add_argument("--seed",type=int,default=42)
    ap.add_argument("--penalizer",type=float,default=0.1)
    args=ap.parse_args()
    df=pd.read_csv(args.data)
    needed={"patient_id","age_years","sex","tumor_grade","survival_days","event"}
    if not needed.issubset(df.columns):
        raise SystemExit(f"Missing required columns: {sorted(needed-set(df.columns))}")
    df=df.dropna(subset=["survival_days","event"]).copy()
    df["survival_days"]=pd.to_numeric(df["survival_days"],errors="coerce")
    df["event"]=pd.to_numeric(df["event"],errors="coerce")
    df=df.dropna(subset=["survival_days","event"])
    df=df[(df["survival_days"]>0)&df["event"].isin([0,1])]
    df["event"]=df["event"].astype(int)
    if len(df)<30 or df["event"].sum()<5:
        raise SystemExit("Insufficient data for this exploratory cross-validation.")
    kf=KFold(n_splits=args.folds,shuffle=True,random_state=args.seed)
    oof=np.full(len(df),np.nan); fold_results=[]
    for fold,(itr,ite) in enumerate(kf.split(df),1):
        train=df.iloc[itr]; test=df.iloc[ite]
        xtr=make_features(train); xte=make_features(test)
        med=xtr.median()
        xtr=xtr.fillna(med); xte=xte.fillna(med)
        scaler=StandardScaler().fit(xtr)
        Xtr=scaler.transform(xtr); Xte=scaler.transform(xte)
        beta=fit_cox(Xtr,train["survival_days"],train["event"],args.penalizer)
        risk=Xte@beta
        oof[ite]=risk
        fold_results.append({"fold":fold,"test_n":int(len(ite)),
          "test_events":int(test["event"].sum()),
          "c_index":concordance_index(test["survival_days"],test["event"],risk)})
    pooled_c=concordance_index(df["survival_days"],df["event"],oof)
    fold_mean=float(np.nanmean([r["c_index"] for r in fold_results]))
    fold_sd=float(np.nanstd([r["c_index"] for r in fold_results],ddof=1))
    X=make_features(df); X=X.fillna(X.median())
    scaler=StandardScaler().fit(X); Xs=scaler.transform(X)
    beta=fit_cox(Xs,df["survival_days"],df["event"],args.penalizer)
    pd.DataFrame({"feature":X.columns,"coefficient_standardized":beta,
                  "hazard_ratio_per_SD":np.exp(beta)}).to_csv(OUT/"clinical_cox_coefficients.csv",index=False)
    pd.DataFrame(fold_results).to_csv(OUT/"cross_validation_folds.csv",index=False)
    metrics={"analysis":"clinical-only exploratory Cox proportional hazards",
      "n_patients":int(len(df)),"n_deaths":int(df["event"].sum()),
      "n_right_censored":int((df["event"]==0).sum()),"features":["age_years","female","high_grade"],
      "cv_folds":int(args.folds),"random_seed":int(args.seed),"penalizer":float(args.penalizer),
      "mean_fold_c_index":fold_mean,"sd_fold_c_index":fold_sd,
      "pooled_out_of_fold_c_index":float(pooled_c),
      "interpretation":"Exploratory estimate on this dataset; small sample and limited predictors. Not clinical validation.",
      "limitation":"No radiomic features were used. Do not claim reproduction of the reference paper's clinical-plus-radiomic model."}
    (OUT/"metrics.json").write_text(json.dumps(metrics,indent=2))
    risk_all=Xs@beta
    plt.figure(figsize=(7,4.5)); plt.hist(risk_all,bins=14)
    plt.xlabel("Relative prognostic score (higher = higher estimated hazard)")
    plt.ylabel("Patients"); plt.title("Clinical-only Cox model risk scores")
    plt.tight_layout(); plt.savefig(OUT/"clinical_risk_scores.png",dpi=180); plt.close()
    groups=oof>=np.nanmedian(oof)
    plt.figure(figsize=(7,4.5))
    for label,mask in [("Lower predicted risk",~groups),("Higher predicted risk",groups)]:
        t=df.loc[mask,"survival_days"].to_numpy(float); e=df.loc[mask,"event"].to_numpy(int)
        unique=np.sort(np.unique(t)); surv=[1.0]; xs=[0.0]; s=1.0
        for tm in unique:
            at_risk=int(np.sum(t>=tm)); deaths=int(np.sum((t==tm)&(e==1)))
            if at_risk and deaths: s*=1-deaths/at_risk
            xs.append(float(tm)); surv.append(s)
        plt.step(xs,surv,where="post",label=label)
    plt.xlabel("Days from diagnosis"); plt.ylabel("Kaplan–Meier survival estimate")
    plt.ylim(0,1.02); plt.title("Survival by out-of-fold clinical risk group")
    plt.legend(); plt.tight_layout(); plt.savefig(OUT/"kaplan_meier_oof_risk_groups.png",dpi=180); plt.close()
    print(json.dumps(metrics,indent=2))
    print(f"Saved model outputs to {OUT}")
if __name__=="__main__": main()
