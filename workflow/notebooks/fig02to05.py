import warnings

import matplotlib.pyplot as plt
import seaborn as sns
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


def plot_group_by_manipulation(df, var_order=None, task_order=None, group_order=None):
    "groupby variable and plot pointplot"
    cmap = sns.color_palette(["#868686", "#D7D7D7"], n_colors=2)
    grid = sns.catplot(
        data=df,
        order=task_order,
        col_order=var_order,
        hue_order=group_order,
        col_wrap=4,
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
    grid.set_xticklabels(rotation=15)
    grid.set_titles(template="{col_name}").set_xlabels("")
    return grid


def _plot_mixedlm(data, **kwargs):
    """plot mixed linear model with random intercept / subject"""

    # estimate coefficients for the 'pain ratings' ~ 'pain signature response'
    variable = data["variable"].iloc[0]
    y = "unp" if "na" in variable else "int"  # unpleasantness or intensity ratings
    mixedlm = smf.mixedlm(
        formula=f"{y} ~ value * group", data=data, groups=data["sub"], missing="drop"
    ).fit()

    # get model predictions
    y_hat = mixedlm.predict(data[["value", "group"]])

    ax = plt.gca()
    ax_kwargs = dict(data=data, x="value", hue="group", palette="binary_r", ax=ax)
    sns.scatterplot(y=y, style="group", alpha=0.75, **ax_kwargs)
    sns.lineplot(y=y_hat, legend=False, lw=3, **ax_kwargs)


def plot_mixedlm(df, var_order=None):
    """plot mixed linear model with random intercept / subject"""
    grid = sns.FacetGrid(df, col="variable", col_order=var_order, col_wrap=4)
    grid.map_dataframe(_plot_mixedlm)
    grid.set_titles(template="{col_name}")
    grid.add_legend()
    grid.set_xlabels("pain signature response").set_ylabels("vas pain rating")
    return grid
