import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from nilearn import image
from nilearn.maskers import NiftiMasker


def replace_path(path, append):
    "replace the extension of a bids path"
    return path.parent / (path.stem.rsplit("_", 4)[0] + append)

def _dot_product(cbf_img, ps_img):
    "compute the dot product of cbf and pain sig images"
    masker = NiftiMasker(image.binarize_img(ps_img), reports=False)
    ps_vect = masker.fit_transform(ps_img)[0]
    cbf_vect = masker.fit_transform(cbf_img)
    return float(np.dot(cbf_vect, ps_vect))

def apply_ps(cbf_img, pain_sigs=["nps", "siips"]):
    "apply each pain sig to the cbf image"
    ps_response = {}
    for ps in pain_sigs:
        ps_img = image.load_img(f"data/pain_sigs/{ps}.nii.gz")
        if ps == "nps":
            cbf_img_resamp = image.resample_to_img(
                cbf_img, ps_img, interpolation="continuous"
            )
            ps_response[ps] = _dot_product(cbf_img_resamp, ps_img)
        else:
            ps_response[ps] = _dot_product(cbf_img, ps_img)
    return ps_response

def main(args):

    ref = pd.read_csv('config/cbf-filtered.csv')
    qc = pd.read_csv(
        replace_path(args.hs, "_desc-quality_control_cbf.csv"),
        usecols=["sub", "ses", "run", "FD", "cbfQEI"],
    )
    df = ref.merge(qc)

    diff = image.math_img("img1 - img2", img1=str(args.hs), img2=str(args.ns))

    if args.smooth_fwhm:
        diff = image.smooth_img(diff, args.smooth_fwhm)

    ps_response = apply_ps(diff, ['nps', 'siips'])
    df.join(pd.DataFrame([ps_response])).to_csv(args.output, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="apply pain signatures to CBF maps")
    parser.add_argument("hs", type=Path, help="input 49C ASLPrep CBF map")
    parser.add_argument("ns", type=Path, help="input 35C ASLPrep CBF map")
    parser.add_argument("output", type=Path, help="output pain signature response")
    parser.add_argument(
        "--smooth_fwhm",
        type=float,
        default=None,
        help="smoothing kernel in mm",
    )

    main(parser.parse_args())
