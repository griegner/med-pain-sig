import warnings

import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import zscore
from statsmodels.stats.multitest import fdrcorrection

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
    data[signature_list] = data[signature_list].dropna().transform(zscore)
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


def plot_group_by_manipulation(df, var_order=None, task_order=None, group_order=None):
    "groupby variable and plot pointplot"
    cmap = sns.color_palette(["#9A9999", "#E9E9E9"], n_colors=2)
    grid = sns.catplot(
        data=df,
        order=task_order,
        col_order=var_order,
        hue_order=group_order,
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
    "plot correlations between pain signature responses and pain ratings"
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


def fit_mixedlm(df, formula=None):
    """mixed linear model with random intercept / subject"""
    mixedlm = smf.mixedlm(formula=formula, data=df, groups=df["sub"], missing="drop").fit()
    return mixedlm.summary().tables[1]


def multicomp_fdr(df, pcol="P>|z|"):
    """the Benjamini-Hochberg procedure"""
    if df.index.get_level_values(1)[0] != "Group Var":
        q_fdr = fdrcorrection(df[pcol].astype("float"))[1]
        df["qFDR"] = q_fdr
    return df
