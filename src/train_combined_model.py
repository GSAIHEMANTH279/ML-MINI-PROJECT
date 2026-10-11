"""Compare clinical-only, radiomics-only, and combined Cox baselines.

Requires data/radiomics_features.csv created by extract_radiomics.py. Feature
selection is performed using training folds only. These exploratory results
are not clinical validation.
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from train_model import fit_cox, concordance_index, make_features

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"outputs"; OUT.mkdir(exist_ok=True)

def prepare_numeric(train_df,test_df,columns):
    tr=train_df[columns].apply(pd.to_numeric,errors="coerce")
    te=test_df[columns].apply(pd.to_numeric,errors="coerce")
    med=tr.median().fillna(0)
    tr=tr.fillna(med).replace([np.inf,-np.inf],np.nan).fillna(med)
    te=te.replace([np.inf,-np.inf],np.nan).fillna(med)
    scaler=StandardScaler().fit(tr)
    return scaler.transform(tr),scaler.transform(te)

def choose_radiomics(train,features,time,event,k=10,penalizer=0.1):
    scores=[]
    for f in features:
        x=pd.to_numeric(train[f],errors="coerce")
        med=x.median()
        if pd.isna(med): continue
        x=x.fillna(med).to_numpy(float)
        if np.nanstd(x)<1e-10: continue
        x=(x-np.mean(x))/(np.std(x)+1e-12)
        try:
            beta=fit_cox(x.reshape(-1,1),time,event,penalizer)
            scores.append((f,abs(float(beta[0]))))
        except Exception:
            continue
    scores.sort(key=lambda z:z[1],reverse=True)
    return [f for f,_ in scores[:min(k,len(scores))]]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--clinical",default=str(ROOT/"data"/"cptac_pda_clinical_survival_clean.csv"))
    ap.add_argument("--radiomics",default=str(ROOT/"data"/"radiomics_features.csv"))
    ap.add_argument("--folds",type=int,default=5)
    ap.add_argument("--seed",type=int,default=42)
    ap.add_argument("--top-k",type=int,default=10)
    ap.add_argument("--penalizer",type=float,default=0.1)
    args=ap.parse_args()
    clinical=pd.read_csv(args.clinical); rad=pd.read_csv(args.radiomics)
    if "patient_id" not in rad.columns: raise SystemExit("Radiomics table needs patient_id column.")
    if rad["patient_id"].duplicated().any():
        raise SystemExit("Radiomics table has duplicate patient_id rows; aggregate/resolve before training.")
    data=clinical.merge(rad,on="patient_id",how="inner",validate="one_to_one")
    if len(data)<30: raise SystemExit(f"Only {len(data)} matched patients; too few for this exploratory 5-fold comparison.")
    radiomic_cols=[c for c in rad.columns if c!="patient_id" and pd.api.types.is_numeric_dtype(pd.to_numeric(rad[c],errors="coerce"))]
    if not radiomic_cols: raise SystemExit("No numeric radiomic feature columns found.")
    clinical_cols=["age_years","female","high_grade"]
    # Derive categorical clinical predictors once using the same baseline definitions.
    clin=make_features(data)
    for c in clinical_cols:
        data[c]=clin[c].to_numpy()
    kf=KFold(n_splits=args.folds,shuffle=True,random_state=args.seed)
    results=[]; oof={m:np.full(len(data),np.nan) for m in ["clinical","radiomics","combined"]}
    for fold,(itr,ite) in enumerate(kf.split(data),1):
        tr=data.iloc[itr].copy(); te=data.iloc[ite].copy()
        selected=choose_radiomics(tr,radiomic_cols,tr["survival_days"].to_numpy(float),tr["event"].to_numpy(int),args.top_k,args.penalizer)
        if not selected: raise SystemExit("No usable radiomic features in a training fold.")
        configs={"clinical":clinical_cols,"radiomics":selected,"combined":clinical_cols+selected}
        for model,cols in configs.items():
            Xtr,Xte=prepare_numeric(tr,te,cols)
            beta=fit_cox(Xtr,tr["survival_days"],tr["event"],args.penalizer)
            risk=Xte@beta; oof[model][ite]=risk
            results.append({"fold":fold,"model":model,"n_train":len(itr),"n_test":len(ite),"selected_radiomics":len(selected),"c_index":concordance_index(te["survival_days"],te["event"],risk)})
    pd.DataFrame(results).to_csv(OUT/"clinical_radiomics_comparison_folds.csv",index=False)
    summary={}
    for model,risk in oof.items():
        foldvals=[r["c_index"] for r in results if r["model"]==model]
        summary[model]={"mean_fold_c_index":float(np.mean(foldvals)),"sd_fold_c_index":float(np.std(foldvals,ddof=1)),"pooled_oof_c_index":float(concordance_index(data["survival_days"],data["event"],risk))}
    report={"n_matched_patients":int(len(data)),"n_radiomic_features_available":int(len(radiomic_cols)),"top_k_selected_per_fold":int(args.top_k),"folds":args.folds,"seed":args.seed,"penalizer":args.penalizer,"results":summary,"important_note":"Radiomic feature selection is performed on training folds only. Results are exploratory and depend on valid masks, extraction settings and matched cases; they are not clinical validation."}
    (OUT/"clinical_radiomics_comparison.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=="__main__": main()
