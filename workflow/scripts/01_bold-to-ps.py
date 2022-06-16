import argparse
from pathlib import Path

import pandas as pd
from nilearn.interfaces import bids

import apply_ps


def get_sub_task(path):
    "get sub and task from a bids path"
    ref = bids.parse_bids_filename(path)
    return ref["sub"], ref["task"]


def main(args):

    sub, task = get_sub_task(args.input)
    ref = (
        pd.read_csv("config/bold-filtered.csv")
        .query(f"sub == '{sub}' and task == '{task}'")
        .reset_index(drop=True)
    )
    ps_response = apply_ps.apply_ps(
        str(args.input),
        ["nps", "siips", "na-gen", "na-therm", "na-mech", "na-sound", "na-vis"],
    )
    ref.join(pd.DataFrame([ps_response])).to_csv(args.output, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="apply pain signatures to pe maps")
    parser.add_argument(
        "input", type=Path, help="input stimulus-evoked parameter estimates"
    )
    parser.add_argument("output", type=Path, help="output pain signature response")

    main(parser.parse_args())
