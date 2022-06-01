import pandas as pd
from sklearn.preprocessing import maxabs_scale


def drop_outliers(df, col):
    "drop outliers if +/- 1.5 IQR"
    q1 = df[col].quantile(0.25)
    q3 = df[col].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    df_drop = df[(df[col] > lower_bound) & (df[col] < upper_bound)]
    print("> dropped outliers:", len(df) - len(df_drop))
    return df_drop


def scale_sig(df):
    cols = df.filter(regex="^sig").columns
    df[cols] = df[cols].transform(maxabs_scale)
    return df


def main():

    (
        pd.read_csv(snakemake.input.path)  # type: ignore
        .assign(group=lambda df: df["sub"].str[0])
        .pipe(drop_outliers, "cbfQEI")
        .pipe(scale_sig)
        .to_csv(snakemake.output[0], index=False)  # type: ignore
    )


if __name__ == "__main__":
    main()
