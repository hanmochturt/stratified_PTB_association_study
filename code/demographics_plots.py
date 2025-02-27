#!/usr/bin/python3
"""
Hannah Takasuka, rotation in the Capra Lab
January 2023
"""

from typing import Dict, List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from tableone import TableOne, load_dataset

import analysis
import constants
import file_format
import stratify_analysis

DEBUG = False  # only use first 50 entries to make the code run faster
PLOT = True
ADJUSTED_ODDS = False

FILENAME = "not_stratified_phenotypes_pvalues.csv"
PHECODES_PER_PATIENT = 'num_phecodes_per_patient'

MERGED_RACE = "Race Hispanic Latino Merged"
SF_CENSUS_2019_RACE = {"White NH": 39.82, "Black NH": 5.23, "Native American NH": 0.28,
                       "Asian NH": 34.57, "Pacific Islander NH": 0.40, "Other NH": 0.41,
                       "Two or more races NH": 4.06, "Hispanic/Latino": 15.24}


def main() -> None:
    """Description"""
    matplotlib.rcParams['font.sans-serif'] = "Arial"  # set default font
    cohort = Cohort()
    # cohort.plot_pie_chart(analysis.RACE)
    # cohort.plot_pie_chart(analysis.PRETERM)
    # plt.figure()
    # cohort.plot_histogram(analysis.GEST_WEEKS)
    # plt.figure()
    # cohort.plot_histogram(analysis.AGE)
    # cohort.plot_early_late_preterm_histograms()
    # cohort.plot_spont_ind_preterm_histograms()
    cohort.plot_spont_ind_box_whisker_violin(analysis.AGE, title="Maternal Age", ylabel="maternal age (years)")
    # cohort.plot_early_late_box_whisker_violin(analysis.AGE, violin=True)
    plot_diagnosis_frequency_histogram(cohort.cohort)
    cohort.plot_spont_ind_box_whisker_violin(PHECODES_PER_PATIENT, outliers=False,
                                             title="Count of Phenotypes per Patient",
                                             ylabel="phenotypes per patient (n)")
    table = cohort.make_table_spont_ind()
    print(table.tabulate(tablefmt="fancy_grid"))
    filename_date_version = file_format.format_filename_with_date("demographics_table.csv", DEBUG)
    table.to_csv(filename_date_version)
    # cohort.plot_early_late_box_whisker_violin(PHECODES_PER_PATIENT, outliers=False)

    plt.show()


class Cohort:

    def __init__(self):
        self.cohort, _, self.phecodes_key = analysis.obtain_data(DEBUG, base_directory=analysis.DIRECTORY)

        print(self.cohort.shape, "***cohort shape")
        self.add_phecodes_per_patient_column()
        self.preterm_cohort = self.cohort[self.cohort[analysis.PRETERM] == True]
        stratified_cohorts = stratify_analysis.stratify_cohort_number_threshold(self.cohort,
                                                                                analysis.GEST_WEEKS,
                                                                                stratify_analysis.EARLY_PRETERM_THRESHOLD_WEEKS)
        self.further_stratified = stratify_analysis.stratify_cohort_number_threshold(
            stratified_cohorts[f"32_and_up"],
            analysis.GEST_WEEKS, 37)

        self.preterm_stratified = {"<32 weeks": stratified_cohorts["less_than_32"],
                                   "32-36 weeks": self.further_stratified["less_than_37"]}

        self.spont_ind_stratified = stratify_analysis.stratify_spont_ind_preterm(self.cohort)

        self.term = self.cohort[self.cohort[analysis.PRETERM] == False]

    def plot_pie_chart(self, column_title: str) -> None:
        dist = analysis.obtain_distribution(column_title, self.cohort)
        fig1, ax1 = plt.subplots()
        ax1.pie(dist.values(), labels=dist.keys(), autopct=make_autopct(dist.values()),
                startangle=90)
        ax1.axis('equal')

    def plot_histogram(self, column_title: str) -> None:
        n, bins, patches = plt.hist(list(self.cohort[column_title]))
        plt.title(column_title)
        plt.xlabel(column_title)
        plt.ylabel("count")

    def plot_early_late_preterm_histograms(self) -> None:
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.suptitle(
            f'{stratify_analysis.CONDITION_LABEL} pre-term pregnancy risk factors identification')

        ax1.hist(list(self.preterm_stratified["<32 weeks"][analysis.AGE]))
        ax1.set_xlabel(f"{analysis.AGE}, <32 weeks")

        ax2.hist(list(self.preterm_stratified["32-36 weeks"][analysis.AGE]))
        ax2.set_xlabel(f"{analysis.AGE}, 32-36 weeks")

    def plot_spont_ind_preterm_histograms(self) -> None:
        cohort_ages = list(self.cohort[analysis.AGE])
        cohort_ages_no_nan = [item for item in cohort_ages if str(item) != 'nan']
        bins = np.histogram(cohort_ages_no_nan, bins=20)[1]
        fig, axs = plt.subplots(3, 1)
        fig.suptitle(
            f'{stratify_analysis.STRATIFY_DATA_VALUES[0]} pre-term pregnancy risk factors identification')

        group_counts = {}
        for index, ax in enumerate(axs):
            if index < 2:
                key = list(self.spont_ind_stratified.keys())[index]
                data = list(self.spont_ind_stratified[key][analysis.AGE])

            else:
                key = "term"
                data = list(self.term[analysis.AGE])
            counts = ax.hist(data, bins)
            group_counts[key] = list(counts[0])
            ax.set_xlabel(f"{analysis.AGE}, {key}")

        plt.figure()
        term_density = [x / len(self.term) for x in list(self.term[analysis.AGE])]
        fig, axs = plt.subplots(2, 1)
        for index, ax in enumerate(axs):
            key = list(self.spont_ind_stratified.keys())[index]
            preterm_stratify_group_density = [float(x) / len(self.term) for x in group_counts[key]]
            preterm_stratify_group_density = [i if i != 0 else 1E-6 for i in preterm_stratify_group_density]
            preterm_over_term_density = [i / j for i, j in zip(term_density, preterm_stratify_group_density)]
            # preterm_over_term_density = [i if i > 1 else -1/i for i in preterm_over_term_density]
            # ax.stairs(preterm_over_term_density, bins)
            middle_of_bin = []
            bins_list = list(bins)
            for i, bin_edge in enumerate(bins_list):
                if i < len(bins_list) - 1:
                    middle_of_bin.append((bin_edge + bins_list[i + 1]) / 2)
            graph = sns.histplot(x=middle_of_bin, bins=bins, weights=preterm_over_term_density, ax=ax,
                                 log_scale=(False, True), element="step", fill=False)
            graph.axhline(1, alpha=0.5)
            ax.set_xlabel(f"{analysis.AGE}, {key} vs. term")
            ax.set_ylabel(f"({key} density this age)/(term density this age)")

    def plot_spont_ind_box_whisker_violin(self, variable: str, ax=None, outliers=True, title=None, ylabel=None) \
            -> None:
        try:
            data = {"indicated": list(self.spont_ind_stratified["indicated"][variable]),
                    "spontaneous": list(self.spont_ind_stratified["spontaneous"][variable]),
                    "term": list(self.term[variable])}
        except KeyError:
            data = {"indicated": list(self.spont_ind_stratified["not_spontaneous"][variable]),
                    "spontaneous": list(self.spont_ind_stratified["spontaneous"][variable]),
                    "term": list(self.term[variable])}
        if ax:
            plot_box_whisker_violin_ax(data, ax, variable, title, ylabel, outliers)
        else:
            plot_box_whisker_violin(data, variable, title, ylabel, outliers)

    def plot_early_late_box_whisker_violin(self, variable: str, outliers=True, violin=False) -> None:
        plt.figure()
        plot_box_whisker_violin({"<32 weeks": list(self.preterm_stratified["<32 weeks"][
                                                       variable]),
                                 "32-36 weeks": list(self.preterm_stratified["32-36 weeks"][
                                                         variable]),
                                 "term": list(self.term[variable])}, variable, violin, outliers)

    def add_phecodes_per_patient_column(self) -> None:
        phecode_sum_term = self.cohort[analysis.get_phecode_columns(self.cohort)].sum(axis=1)
        self.cohort[PHECODES_PER_PATIENT] = phecode_sum_term

    def make_table_spont_ind(self) -> TableOne:
        df_combined = stratify_analysis.move_dict_group_key_to_column(
            {**self.spont_ind_stratified, **{"term": self.term}})
        print(analysis.obtain_distribution(stratify_analysis.STRATIFY_GROUP, df_combined))
        df_combined = add_columns_all_preterm_and_overall(df_combined)
        print(analysis.obtain_distribution(stratify_analysis.STRATIFY_GROUP, df_combined))

        df_combined_edu = analysis.add_maternal_edu_simplified_column(df_combined)
        median_min_max = [analysis.AGE, PHECODES_PER_PATIENT]
        mean_sd_titles = []
        for column_title in median_min_max:
            mean_sd_title = f"{column_title}_sd"
            df_combined_edu[mean_sd_title] = df_combined_edu[column_title]
            mean_sd_titles.append(mean_sd_title)
        table = TableOne(df_combined_edu, categorical=[analysis.RACE, analysis.INS_PRIVATE, analysis.EDU_SIMP],
                         columns=[analysis.RACE, analysis.AGE, analysis.INS_PRIVATE, analysis.EDU_SIMP,
                                  PHECODES_PER_PATIENT] + mean_sd_titles,
                         groupby=[stratify_analysis.STRATIFY_GROUP], overall=False, nonnormal=median_min_max,
                         min_max=median_min_max)
        return table


def add_columns_all_preterm_and_overall(df_spont_ind: pd.DataFrame) -> pd.DataFrame:
    df = df_spont_ind.copy()
    additional_df_columns = {"all preterm": df_spont_ind[df_spont_ind[analysis.PRETERM] == True],
                             "overall": df_spont_ind}
    for title, sub_df in additional_df_columns.items():
        sub_df[stratify_analysis.STRATIFY_GROUP] = title
        #df = df.append(sub_df, ignore_index=True)
        df = pd.concat([df, sub_df], ignore_index=True)
    return df


def plot_box_whisker_violin(data: Dict, variable="", title=None, ylabel=None, outliers=False) -> None:
    plt.figure()
    ax = plt.axes()
    plot_box_whisker_violin_ax(data, ax, variable, title, ylabel, outliers)


def plot_box_whisker_violin_ax(data: Dict, ax: matplotlib.axes.Axes, variable="", title=None, ylabel=None,
                               outliers=False) -> None:
    data_no_nan = analysis.remove_na_from_dict(data)

    if not title:
        title = variable
    if not ylabel:
        ylabel = variable
    # vp = ax.violinplot(data_no_nan.values(), widths=0.5, showextrema=outliers)
    bp = ax.boxplot(data_no_nan.values(), widths=0.5, patch_artist=True, showfliers=outliers)  # sym=outlier_symbol)

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    # Label x-axis ticks
    xticklabels = [label.title() for label in list(data_no_nan.keys())]
    ax.set_xticklabels(xticklabels)
    xticks = [x + 1 for x in range(len(xticklabels))]
    ax.set_xticks(xticks)

    for index, patch in enumerate(bp['boxes']):
        patch.set(facecolor=constants.COLORS[xticklabels[index]])
    for element in ['medians', ]:  # 'whiskers', 'means',  'caps']:
        plt.setp(bp[element], color='black')

    add_significance_box_whisker(ax, data)


def add_significance_box_whisker(ax: matplotlib.axes.Axes, data: Dict) -> None:
    bottom, top = ax.get_ylim()
    y_range = top - bottom
    data_no_nan = analysis.remove_na_from_dict(data)
    significant_combinations = calculate_significance(list(data_no_nan.values()))
    for i, significant_combination in enumerate(significant_combinations):
        # Columns corresponding to the datasets of interest
        x1 = significant_combination[0][0]
        x2 = significant_combination[0][1]
        # What level is this bar among the bars above the plot?
        level = len(significant_combinations) - i
        # Plot the bar
        bar_height = (y_range * 0.07 * level) + top
        bar_tips = bar_height - (y_range * 0.02)
        ax.plot(
            [x1, x1, x2, x2],
            [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
        )
        # Significance level
        p = significant_combination[1]
        if p < 0.001:
            sig_symbol = f'***, p={analysis.round_sig_figs(p)}'
        elif p < 0.01:
            sig_symbol = f'**, p={analysis.round_sig_figs(p)}'
        elif p < 0.05:
            sig_symbol = f'*, p={analysis.round_sig_figs(p)}'
        else:
            sig_symbol = f'ns, p={analysis.round_sig_figs(p)}'
        text_height = bar_height + (y_range * 0.01)
        ax.text((x1 + x2) * 0.5, text_height, sig_symbol, ha='center', va='bottom', c='k')


def calculate_significance(data: list) -> List:
    # Initialise a list of combinations of groups that are significantly different
    significant_combinations = []
    # Check from the outside pairs of boxes inwards
    ls = list(range(1, len(data) + 1))
    combinations = [(ls[x], ls[x + y]) for y in reversed(ls) for x in range((len(ls) - y))]
    for combination in combinations:
        data1 = data[combination[0] - 1]
        data2 = data[combination[1] - 1]
        # Significance
        U, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        # if p < 0.05:
        significant_combinations.append([combination, p])
    return significant_combinations


def plot_diagnosis_frequency_histogram_x(cohort: pd.DataFrame) -> None:
    plt.figure()
    phecode_sum = measure_diagnosis_count(cohort)
    phe_list = list(phecode_sum)
    phe_list.sort()
    n, bins, patches = plt.hist(list(phecode_sum))
    plt.title(f"Diagnosis Frequency Distribution")
    plt.xlabel("# of patients with diagnosis")
    plt.ylabel("# of diagnoses")


def plot_diagnosis_frequency_histogram(cohort: pd.DataFrame, ax=None) -> None:
    if ax:
        ax.set_title("Diagnosis Frequency Distribution")
        ax.set_xlabel("# of patients with diagnosis")
        ax.set_ylabel("# of diagnoses")
    else:
        plt.figure()
        plt.title("Diagnosis Frequency Distribution")
        plt.xlabel("# of patients with diagnosis")
        plt.ylabel("# of diagnoses")

    phecode_sum = measure_diagnosis_count(cohort)
    phe_list = list(phecode_sum)
    phe_list.sort()
    splot = sns.histplot(list(phecode_sum), log_scale=(False, True), binwidth=100, ax=ax)



def measure_diagnosis_count(cohort: pd.DataFrame) -> pd.DataFrame:
    phecode_sum = cohort[analysis.get_phecode_columns(cohort)].sum(axis=0)
    return phecode_sum


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return '{p:.0f}%  ({v:d})'.format(p=pct, v=val)

    return my_autopct


if __name__ == "__main__":
    main()
