import pandas as pd
from sklearn.preprocessing import maxabs_scale


def scale_sig(df):
    cols = df.filter(regex="^sig").columns
    df[cols] = df[cols].transform(maxabs_scale)
    return df


def main():

    (
        pd.read_csv(snakemake.input.path)  # type: ignore
        .assign(group=lambda df: df["sub"].str[0])
        .query("FD < 0.5")
        .pipe(scale_sig)
        .to_csv(snakemake.output[0], index=False)  # type: ignore
    )


if __name__ == "__main__":
    main()
