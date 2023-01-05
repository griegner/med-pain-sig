import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from nilearn.glm.first_level import FirstLevelModel
from nilearn.interfaces import bids, fmriprep


class GLMInputs:
    "mask, events, and confounds attributes of the preproc input"

    def __init__(self, preproc):
        self.preproc = str(preproc)
        self.ref = bids.parse_bids_filename(preproc)
        self.mask = str(preproc.parent / f"{preproc.name[:-20]}-brain_mask.nii.gz")
        self.events = bids.get_bids_files(
            preproc.parents[4] / "rawdata",
            file_tag="events",
            file_type="tsv",
            modality_folder="func",
            sub_label=self.ref["sub"],
            filters=[("task", self.ref["task"])],
        )[0]
        self.confounds, self.sample_mask = fmriprep.load_confounds(
            str(preproc),
            strategy=("motion", "high_pass", "compcor", "scrub"),
            motion="basic",
            compcor="anat_combined",
            n_compcor=5,
            scrub=0,
            fd_threshold=0.5,
            std_dvars_threshold=3,
        )
        self.confounds = self.add_scrub()

    def add_scrub(self):
        "add scrubbing columns to confounds"
        if self.sample_mask is None:
            return self.confounds
        else:
            n_volumes = len(self.confounds)
            sample_mask_index = np.setdiff1d(np.arange(n_volumes), self.sample_mask)
            sample_mask_one_hot = np.zeros((n_volumes, len(sample_mask_index)))
            sample_mask_one_hot[
                sample_mask_index, np.arange(len(sample_mask_index))
            ] = 1
            return pd.concat(
                [self.confounds, pd.DataFrame(sample_mask_one_hot)], axis=1
            )


def main(args):

    inputs = GLMInputs(args.input)

    glm = FirstLevelModel(
        t_r=2.0,
        hrf_model="glover + derivative",
        drift_model=None,
        mask_img=inputs.mask,
        standardize=True,
        signal_scaling=False,
        smoothing_fwhm=args.smooth_fwhm,
        n_jobs=args.n_jobs,
    ).fit(run_imgs=inputs.preproc, events=inputs.events, confounds=inputs.confounds)
    pe = glm.compute_contrast("stim", output_type="effect_size")
    pe.to_filename(args.output)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="apply pain signatures to BOLD maps")
    parser.add_argument("input", type=Path, help="input fMRIprep preprocessed data")
    parser.add_argument(
        "output", type=Path, help="output stimulus-evoked parameter estimates"
    )
    parser.add_argument(
        "--smooth_fwhm",
        type=float,
        default=None,
        help="smoothing kernel in mm",
    )
    parser.add_argument("--n_jobs", type=int, default=1, help="n_jobs")

    main(parser.parse_args())
