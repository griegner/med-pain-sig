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
        replace_path(args.input, "_desc-quality_control_cbf.csv"),
        usecols=["sub", "ses", "run", "FD", "cbfQEI"],
    )

    img = image.load_img(str(args.input))

    if args.smooth_fwhm:
        img = image.smooth_img(img, args.smooth_fwhm)

    ps_response = apply_ps.apply_ps(img, ["nps", "siips", "pines"])

    df = ref.merge(qc)
    df.join(pd.DataFrame([ps_response])).to_csv(args.output, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="apply pain signatures to CBF maps")
    parser.add_argument("input", type=Path, help="input ASLPrep CBF map")
    parser.add_argument("output", type=Path, help="output pain signature response")
    parser.add_argument(
        "--smooth_fwhm",
        type=float,
        default=None,
        help="smoothing kernel in mm",
    )

    main(parser.parse_args())
