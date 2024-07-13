import matplotlib.pyplot as plt
import seaborn as sns
from pymer4.models import Lmer
from scipy.stats import zscore
from statsmodels.stats.multitest import fdrcorrection

plt.rcParams["font.size"] = 18


def scale(data, signature_list):
    "scale signature responses to zero mean and unit variance"
    data[signature_list] = data[signature_list].dropna().transform(zscore)
    return data


def plot_group_by_manipulation(df, var_order=None, task_order=None, group_order=None, ylim=None):
    "groupby variable and plot pointplot"
    cmap = sns.color_palette(["#969191", "#ECE6D0"], n_colors=2)
    grid = sns.catplot(
        data=df,
        order=group_order,
        col_order=var_order,
        hue_order=task_order,
        kind="bar",
        x="group",
        y="value",
        col="variable",
        hue="task",
        palette=cmap,
        errorbar="se",
        capsize=0.15,
        aspect=2,
        edgecolor="k",
        linewidth=0.2,
        sharey=False,
    )
    grid.set(ylim=ylim)
    grid.set_titles(template="{col_name}").set_xlabels("")
    [ax.axhline(y=0, c="k", lw=1) for ax in grid.axes.flat]
    return grid


def _get_sig_labels(pval):
    """categorize p-values by significance labels"""
    if pval <= 0.001:
        return "***"
    elif pval <= 0.01:
        return "**"
    elif pval <= 0.05:
        return "*"
    else:
        return ""


def anova_table(df, formula=None, tasks=None, groups=None):
    """ANOVA table from a fitted lme4 model"""
    mixedlm = Lmer(formula=formula, data=df)
    mixedlm.fit(factors={"task": tasks, "group": groups}, summarize=False, REML=False)
    return mixedlm.anova()


def paired_ttests(df, formula=None, tasks=None, groups=None):
    """paired t-tests between marginal means, p-value correction by fdr"""
    mixedlm = Lmer(formula=formula, data=df)
    mixedlm.fit(
        factors={"task": tasks, "group": groups, "gender": ["m", "f"]}, summarize=False, REML=False
    )
    marginals, ttests = mixedlm.post_hoc(
        marginal_vars=["task"], grouping_vars=["group"], p_adjust=None
    )
    ttests = ttests.query("~Contrast.str.contains('baseline')")
    ttests.insert(ttests.shape[1] - 1, "P-fdr", fdrcorrection(ttests["P-val"])[1])
    ttests["Sig"] = ttests["P-fdr"].apply(_get_sig_labels)
    return ttests


def unpaired_ttests(df, formula=None, groups=None):
    """unpaired t-tests between marginal differences, p-value correction by fdr"""
    mixedlm_diff = Lmer(formula=formula, data=df)
    mixedlm_diff.fit(factors={"group": groups, "gender": ["m", "f"]}, summarize=False, REML=False)
    marginals, ttests = mixedlm_diff.post_hoc(marginal_vars=["group"], p_adjust=None)
    ttests.insert(8, "P-fdr", fdrcorrection(ttests["P-val"])[1])
    ttests["Sig"] = ttests["P-fdr"].apply(_get_sig_labels)
    return ttests
