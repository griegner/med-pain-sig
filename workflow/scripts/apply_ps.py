import numpy as np
from nilearn import image
from nilearn.maskers import NiftiMasker


def _dot_product(img, ps_img):
    "compute the dot product of input and pain signature images"
    masker = NiftiMasker(image.binarize_img(ps_img), reports=False)
    ps_vect = masker.fit_transform(ps_img)[0]
    img_vect = masker.fit_transform(img)
    return float(np.dot(img_vect, ps_vect))


def apply_ps(img, pain_sigs=["nps", "siips"]):
    "apply each pain signature to the input image"
    ps_response = {}
    for ps in pain_sigs:
        ps_img = image.load_img(f"data/pain_sigs/{ps}.nii.gz")
        if ps == "nps":
            img_resamp = image.resample_to_img(img, ps_img, interpolation="continuous")
            ps_response[ps] = _dot_product(img_resamp, ps_img)
        else:
            ps_response[ps] = _dot_product(img, ps_img)
    return ps_response
