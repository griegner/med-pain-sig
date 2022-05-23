import argparse
from pathlib import Path

import pandas as pd
from nilearn import image

import apply_ps


def replace_path(path, append):
    "replace the extension of a bids path"
    return path.parent / (path.stem.rsplit("_", 4)[0] + append)


def main(args):

    ref = pd.read_csv("config/cbf-filtered.csv")
    qc = pd.read_csv(
        replace_path(args.hs, "_desc-quality_control_cbf.csv"),
        usecols=["sub", "ses", "run", "FD", "cbfQEI"],
    )
    df = ref.merge(qc)

    diff = image.math_img("img1 - img2", img1=str(args.hs), img2=str(args.ns))

    if args.smooth_fwhm:
        diff = image.smooth_img(diff, args.smooth_fwhm)

    ps_response = apply_ps.apply_ps(diff, ["nps", "siips"])
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
