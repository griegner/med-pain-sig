import warnings

import matplotlib.pyplot as plt
import pingouin as pg
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import zscore

warnings.simplefilter("ignore", category=FutureWarning)


def drop_outliers(df, col):
    "drop outliers if +/- 1.5 IQR"
    q1 = df[col].quantile(0.25)
    q3 = df[col].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    return df[df[col] > lower_bound]


def scale(data, signature_list):
    "scale signature responses to zero mean and unit variance"
    data[signature_list] = data[signature_list].transform(zscore)
    return data


def get_residual(df):
    "return residual intensity controlling for nps response"
    y = df["nps"]
    X = df["int"]
    X = sm.add_constant(X)
    df["int_residual"] = sm.OLS(y, X).fit().resid
    return df


def get_corrwith_pain(df):
    """pearson correlations: na-unp, nps-int, siips-int_residual"""
    signature = df["variable"].unique()[0]
    if "na" in signature:
        return df[["value"]].corrwith(df["unp"])
    elif "nps" in signature:
        return df[["value"]].corrwith(df["int"])
    elif "siips" in signature:
        return df[["value"]].corrwith(df["int_residual"])


def plot_group_by_manipulation(df, order=None, col_order=None):
    "groupby variable and plot pointplot"
    cmap = sns.color_palette(["#E9E9E9", "#9A9999"], n_colors=2)
    grid = sns.catplot(
        data=df,
        order=order,
        col_order=col_order,
        kind="bar",
        x="task",
        y="value",
        col="variable",
        hue="group",
        palette=cmap,
        errorbar="se",
        aspect=1,
        height=3,
        edgecolor="k",
        linewidth=0.2,
    )
    (grid.set_xticklabels(rotation=15).set_titles(template="{col_name}").set_xlabels(""))
    return grid


def plot_correlations(df, order=None):
    "plot correlations between brain pattern responses and pain ratings"
    cmap = sns.color_palette(["#D62F26", "#EF633E"] + 5 * ["#4067AD"], n_colors=7)
    ax = sns.pointplot(
        data=df,
        x="variable",
        y="value",
        errorbar="se",
        palette=cmap,
        order=order,
        join=False,
    )
    ax.axhline(0, linestyle="--", lw=0.5, color="k")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15)
    ax.set_ylabel("correlation (r) with pain")
    ax.set_xlabel("")
    ax.set_yticks([-0.1, 0, 0.1, 0.2, 0.3, 0.4])
    plt.show()
    return ax.get_figure()


def _multicomp_fdr(df, pcol="p-unc"):
    """the Benjamini-Hochberg procedure"""
    p_unc = df[pcol]
    # returns FDR-corrected p-values
    q_fdr = pg.multicomp(p_unc, method="fdr_bh")[1]
    idx = df.columns.get_loc(pcol)
    df.insert(idx + 1, "q-fdr", q_fdr)
    return df


def mixed_anova(df, alpha=0.05):
    "group x manipulation anova"

    kwargs = dict(dv="value", within="task", between="group", subject="sub")

    mixed_anova_df = (
        df.groupby("variable")
        .apply(lambda df: df.mixed_anova(**kwargs))
        .groupby("Source")
        .apply(_multicomp_fdr)
    )

    sig_interaction = mixed_anova_df.query(
        "Source == 'Interaction' and `q-fdr` < @alpha"
    ).index.get_level_values(0)

    pairwise_ttests_df = (
        df.query("variable in @sig_interaction")
        .groupby("variable")
        .apply(
            lambda df: df.pairwise_tests(
                within_first=False, nan_policy="pairwise", effsize="cohen", **kwargs
            )
        )
        .query("Contrast == 'group * task'")
        .pipe(_multicomp_fdr)
    )

    return mixed_anova_df, pairwise_ttests_df


def one_sample_ttest(df):
    "one sample ttest"
    return (
        df.groupby(["group", "variable"])
        .apply(lambda df: pg.ttest(df["value"], 0))
        .pipe(_multicomp_fdr, "p-val")
    )
