"""Extract PyRadiomics features from CT volumes and matching tumor masks.

Prepare data/radiomics_manifest.csv with columns:
patient_id,image_path,mask_path
Paths must point to a medical image volume and a matching segmentation mask
(e.g. NIfTI files). This script does not download DICOMs or create tumor masks.
"""
from pathlib import Path
import argparse
import pandas as pd

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",default="data/radiomics_manifest.csv")
    ap.add_argument("--output",default="data/radiomics_features.csv")
    ap.add_argument("--bin-width",type=float,default=25.0)
    args=ap.parse_args()
    try:
        from radiomics import featureextractor
    except ImportError as exc:
        raise SystemExit("Install PyRadiomics and SimpleITK first: pip install pyradiomics SimpleITK") from exc
    manifest=pd.read_csv(args.manifest)
    required={"patient_id","image_path","mask_path"}
    if not required.issubset(manifest.columns):
        raise SystemExit(f"Manifest must include columns: {sorted(required)}")
    extractor=featureextractor.RadiomicsFeatureExtractor(binWidth=args.bin_width)
    extractor.enableAllImageTypes()
    rows=[]; errors=[]
    for _,row in manifest.iterrows():
        pid=str(row["patient_id"])
        image_path=Path(str(row["image_path"]))
        mask_path=Path(str(row["mask_path"]))
        if not image_path.is_file() or not mask_path.is_file():
            errors.append({"patient_id":pid,"error":"image or mask path does not exist"})
            continue
        try:
            values=extractor.execute(str(image_path),str(mask_path))
            features={k:v for k,v in values.items() if not str(k).startswith("diagnostics_")}
            rows.append({"patient_id":pid,**features})
        except Exception as exc:
            errors.append({"patient_id":pid,"error":str(exc)})
    if not rows:
        raise SystemExit("No radiomics rows extracted. Check image/mask paths and label value.")
    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out,index=False)
    if errors:
        err_path=out.with_name("radiomics_extraction_errors.csv")
        pd.DataFrame(errors).to_csv(err_path,index=False)
        print(f"Rows with errors: {len(errors)}; details: {err_path}")
    print(f"Saved {len(rows)} patients and {len(pd.DataFrame(rows).columns)-1} feature columns to {out}")
    print("Check the mask label value, image orientation/spacing, and patient IDs before survival modelling.")

if __name__=="__main__":
    main()
