import matplotlib.pyplot as plt
import pingouin as pg
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import zscore


def drop_outliers(df, col):
    "drop outliers if +/- 1.5 IQR"
    q1 = df[col].quantile(0.25)
    q3 = df[col].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    return df[df[col] > lower_bound]


def scale(data, signature_list):
    "scale signature reponses to zero mean and unit variance"
    data[signature_list] = data[signature_list].transform(zscore)
    return data


def get_residual(df):
    "return residual intensity controlling for nps response"
    y = df["nps"]
    X = df["int"]
    X = sm.add_constant(X)
    return sm.OLS(y, X).fit().resid


def get_corrwith_pain(df):
    signature = df["variable"].unique()[0]
    if "na" in signature:
        return df[["value"]].corrwith(df["unp"])
    elif "nps" in signature:
        return df[["value"]].corrwith(df["int"])
    elif "siips" in signature:
        return df[["value"]].corrwith(df["int_residual"])


def plot_group_by_manipulation(df, order=None, col_order=None, join=True):
    "groupby variable and plot pointplot"
    grid = sns.catplot(
        data=df,
        order=order,
        col_order=col_order,
        join=join,
        kind="point",
        x="task",
        y="value",
        col="variable",
        hue="group",
        palette="Blues",
        dodge=True,
        errorbar="se",
        aspect=1,
        height=3,
    )
    (
        grid.set_xticklabels(rotation=15)
        .set_titles(template="{col_name}")
        .set_xlabels("")
    )
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


def mixed_anova(df, alpha=0.05):
    "group x manipulation anova"

    kwargs = dict(dv="value", within="task", between="group", subject="sub")

    mixed_anova_df = df.groupby("variable").apply(lambda df: df.mixed_anova(**kwargs))

    sig_interaction = (
        mixed_anova_df.query("Source == 'Interaction'")
        .pipe(
            lambda df: df[
                pg.multicomp(df["p-unc"].values, alpha=alpha, method="fdr_bh")[0]
            ]
        )  # multiple comparisons
        .index.get_level_values(0)
    )

    pairwise_ttests_df = (
        df.query("variable in @sig_interaction")
        .groupby("variable")
        .apply(
            lambda df: df.pairwise_tests(
                within_first=False, nan_policy="pairwise", **kwargs
            )
        )
        .query("Contrast == 'group * task'")
    )

    return mixed_anova_df, pairwise_ttests_df


def one_sample_ttest(df):
    "one sample ttest"
    return df.groupby("variable").apply(lambda df: pg.ttest(df["value"], 0))
