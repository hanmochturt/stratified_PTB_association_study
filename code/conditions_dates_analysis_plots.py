#!/usr/bin/python3
"""
Hannah Takasuka, rotation in the Capra Lab
January 2023
"""
import datetime
import time
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import numpy as np
import pandas as pd
from tableone import TableOne

import analysis
import constants
import demographics_plots as dp
import file_format
import stratify_analysis

DEBUG = False  # only use first 50 entries to make the code run faster
COND_DATES = "conditions_raw"

#df labels: (can be edited)
FIRST_VISIT = "first visit, years before conception"
ALL_VISIT_DATES = "diagnoses days before conception"


def main() -> None:
    matplotlib.rcParams['font.sans-serif'] = "Arial"  # set default font
    cohort = get_cohort_patients_conds_dates()

    #plot_dates_box_whisker_stratified(cohort, "total_diagnoses")

    plot_dates_box_whisker_stratified(cohort, FIRST_VISIT,
                                    title="Time in EHR before Conception",
                                    ylabel="years before conception")

    plot_diagnoses_time_before_conception_stratified(cohort)

    #table = make_table_spont_ind(cohort)
    #print(make_table_spont_ind(cohort).tabulate(tablefmt="fancy_grid"))
    print(make_table_spont_ind_visit_dates(cohort).tabulate(tablefmt="fancy_grid"))
    plt.show()


def get_cohort_patients_conds_dates() -> pd.DataFrame:
    cohort, _, phecodes_key = analysis.obtain_data(DEBUG, analysis.DIRECTORY)
    conditions_dates = get_conds_dates()
    cohort_conditions_filtered = filter_for_overlapping_infant_ids(conditions_dates, cohort)
    person_diagnosis_dates = get_visits_with_diagnosis_per_patient(cohort_conditions_filtered)
    cohort[ALL_VISIT_DATES] = person_diagnosis_dates
    cohort[FIRST_VISIT] = cohort[ALL_VISIT_DATES].map(lambda x: max(x).days / 365.25)
    cohort_num_diagnoses = add_total_diagnoses_column(conditions_dates, cohort)
    cohort_num_diagnoses = stratify_analysis.add_spont_ind_term_column(cohort_num_diagnoses)
    return cohort_num_diagnoses


def get_conds_dates() -> pd.DataFrame:
    directories = analysis.set_directories(analysis.DIRECTORY)
    filename = file_format.find_recent_datetime_file(COND_DATES, directories[analysis.DATA_FOLDER])
    print(filename)
    print(f"raw conditions file uploading. This takes ~5-13 minutes. Started at {datetime.datetime.now().time()}")
    start_time = time.time()
    conditions_dates = pd.read_csv(filename, encoding="cp1252")
    print(f"raw conditions file uploaded. {(time.time()-start_time)//60} minutes elapsed.")
    return conditions_dates


def filter_for_overlapping_infant_ids(conditions_dates: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    cohort_person_ids_inf = list(cohort.index)
    cohort_conditions = conditions_dates.loc[conditions_dates[analysis.ID_INF].isin(cohort_person_ids_inf)]
    return cohort_conditions


def get_visits_with_diagnosis_per_patient(df: pd.DataFrame) -> Dict:
    df['condition_days_before_conception'] = pd.to_datetime(df["preg_start"]) - pd.to_datetime(
        df['condition_start_date'])
    df_no_duplicates = df.drop_duplicates(subset=[analysis.ID_INF, 'condition_days_before_conception'], keep='first')
    df_no_duplicates = df_no_duplicates.set_index(analysis.ID_INF)
    return df_no_duplicates.groupby(level=0)['condition_days_before_conception'].apply(list)


def reformat_visits_per_patient_to_years_list(cohort: pd.DataFrame) -> List:
    person_diagnosis_dates = cohort["diagnoses days before conception"].to_dict()
    diagnoses_days_before_pregnancy_2d = list(person_diagnosis_dates.values())
    diagnoses_days_before_pregnancy_dt = [item for sublist in diagnoses_days_before_pregnancy_2d for item in sublist]
    diagnoses_days_before_pregnancy = [i.days for i in diagnoses_days_before_pregnancy_dt]
    diagnoses_years_before_pregnancy = [x / 365.25 for x in diagnoses_days_before_pregnancy]
    return diagnoses_years_before_pregnancy


#def merge_diganoses_dates_and_cohort(person_diagnosis_dates: Dict, cohort:pd.DataFrame) -> pd.DataFrame:


def plot_years_histogram(years: List) -> None:
    fig, ax = plt.subplots()
    ax.hist(years)
    ax.set_xlabel("Years Before Conception")
    ax.set_ylabel("Number of visits with a diagnosis")
    ax.invert_xaxis()


def add_total_diagnoses_column(conditions_dates: pd.DataFrame, cohort_no_num_diag: pd.DataFrame) -> pd.DataFrame:
    personid_num_diagnoses = analysis.obtain_distribution("person_id_inf", conditions_dates)
    df_num_diagnoses = pd.Series(personid_num_diagnoses, name="total_diagnoses")
    cohort_with_num_diagnoses = pd.concat(
        [df_num_diagnoses, cohort_no_num_diag], axis=1, join='inner')
    cohort_with_num_diagnoses = analysis.define_preterm(cohort_with_num_diagnoses)
    return cohort_with_num_diagnoses


def plot_dates_box_whisker_stratified(cohort, variable, outliers=False, violin=False, title=None, ylabel=None) -> None:
    stratified_cohort = stratify_analysis.stratify_spont_ind_preterm(cohort, term_separate=True)
    dp.plot_box_whisker_violin({"indicated": list(stratified_cohort["indicated"][
                                                      variable]),
                                "spontaneous": list(stratified_cohort["spontaneous"][
                                                        variable]),
                                "term": list(stratified_cohort["term"][variable])}, variable, violin, title, ylabel,
                               outliers)


def plot_diagnoses_time_before_conception_stratified_hist(cohort_with_diagnosis_dates: pd.DataFrame, log=True, ax=None) -> None:
    if ax:
        ax.set_title("Diagnosis Time Distribution")
        ax.set_xlabel("Years before Conception")
        ax.set_ylabel("Number of visits with a diagnosis (density)")
    else:
        plt.figure()
        plt.xlabel("Years before Conception")
        plt.ylabel("Number of visits with a diagnosis (density)")
    group_titles = analysis.obtain_distribution(stratify_analysis.STRATIFY_GROUP, cohort_with_diagnosis_dates)
    '''group_visit_days_per_patient_dict = {**dict.fromkeys((group_titles.keys()), False)}
    for group in group_titles.keys():
        group_visit_days_per_patient_dict[group] = cohort_with_diagnosis_dates.loc[
            cohort_with_diagnosis_dates[stratify_analysis.STRATIFY_GROUP] == group]
        group_visit_days_per_patient_dict[group] = group_visit_days_per_patient_dict[group][
            "diagnoses days before conception"]'''

    group_visit_days_dict = {**dict.fromkeys((group_titles.keys()), False)}
    all_days = []

    for group in group_titles.keys():
        group_visit_days_dict[group] = reformat_visits_per_patient_to_years_list(
            cohort_with_diagnosis_dates.loc[
                cohort_with_diagnosis_dates[stratify_analysis.STRATIFY_GROUP] == group])
        all_days += group_visit_days_dict[group]

    bins = np.histogram(all_days, bins=24)[1]
    for group in group_titles.keys():
        if not ax:
            plt.hist(group_visit_days_dict[group], bins=bins, histtype='bar', density=True,
                    label=group, alpha=0.75, log=log, color=constants.COLORS[group.capitalize()])
            plt.legend()
            lower_x, upper_x = plt.xlim()
            plt.xlim(upper_x, lower_x)

    if ax:
        colors = []
        for group in group_visit_days_dict.keys():
            colors.append(constants.COLORS[group.capitalize()])
        ax.hist(group_visit_days_dict.values(), bins=bins, histtype='bar', density=True,
                label=list(group_visit_days_dict.keys()), log=log, color=colors)
        lower_x, upper_x = ax.get_xlim()
        ax.set_xlim(upper_x, lower_x)
        ax.yaxis.tick_right()
    ax2 = plt.axes([0, 0, 1, 1])
    ax2.yaxis.tick_right()
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax, [0.02, 0.2, 0.61, 0.7])
    ax2.set_axes_locator(ip)
    #axins = inset_axes(ax, ax2, loc1=2, loc2=4, bbox_to_anchor=(.08, 0.35))#, bbox_transform=axfft.figure.transFigure)

    ax2.hist(group_visit_days_dict.values(), bins=bins, histtype='bar', density=True,
                label=list(group_visit_days_dict.keys()), log=log, color=colors)
    ax2.legend()
    patch, pp1,pp2 = mark_inset(ax, ax2, loc1=1, loc2=1, linewidth=0.75, fc="None", ec='k', alpha=0.4, clip_on=True, zorder=1)

    pp1.loc1 = 3
    pp1.loc2 = 4
    pp2.loc1 = 4
    pp2.loc2 = 3

    x1, x2, y1, y2 = 22.2, 8, 0, 0.003  # specify the limits
    ax2.set_xlim(x1, x2)  # apply the x-limits
    ax2.set_ylim(y1, y2)  # apply the y-limits
    ax2.set_xticks([])

    #mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")





def plot_diagnoses_time_before_conception_stratified(cohort_with_diagnosis_dates: pd.DataFrame) -> None:
    plt.figure()
    group_titles = analysis.obtain_distribution(stratify_analysis.STRATIFY_GROUP, cohort_with_diagnosis_dates)
    group_visit_days_dict = {**dict.fromkeys((group_titles.keys()), False)}
    all_days = []

    for group in group_titles.keys():
        group_visit_days_dict[group] = reformat_visits_per_patient_to_years_list(
            cohort_with_diagnosis_dates.loc[
                cohort_with_diagnosis_dates[stratify_analysis.STRATIFY_GROUP] == group])
        all_days += group_visit_days_dict[group]

    dp.plot_box_whisker_violin({"indicated": list(group_visit_days_dict["indicated"]),
                        "spontaneous": list(group_visit_days_dict["spontaneous"]),
                        "term": list(group_visit_days_dict["term"])}, ylabel="Years before Conception", title="Time of Visits with a Diagnosis")


def make_table_spont_ind(df_spont_ind: pd.DataFrame) -> TableOne:
    df = dp.add_columns_all_preterm_and_overall(df_spont_ind).copy()
    df_combined_edu = analysis.add_maternal_edu_simplified_column(df)
    median_min_max = [FIRST_VISIT] #TODO add all visits
    mean_sd_titles = []
    for column_title in median_min_max:
        mean_sd_title = f"{column_title}_sd"
        df_combined_edu[mean_sd_title] = df_combined_edu[column_title]
        mean_sd_titles.append(mean_sd_title)
    table = TableOne(df_combined_edu,
                     columns=[FIRST_VISIT]+mean_sd_titles,
                     groupby=[stratify_analysis.STRATIFY_GROUP], overall=False, nonnormal=median_min_max,
                     min_max=median_min_max)
    return table


def make_table_spont_ind_visit_dates(df_spont_ind: pd.DataFrame) -> TableOne:
    df = dp.add_columns_all_preterm_and_overall(df_spont_ind)
    df = df.explode(ALL_VISIT_DATES, ignore_index=True)
    df[ALL_VISIT_DATES] = df[ALL_VISIT_DATES].map(lambda x: x.days/365.25)
    column_title = ALL_VISIT_DATES
    mean_sd_title = f"{column_title}_sd"
    df[mean_sd_title] = df[column_title]
    table = TableOne(df,
                     columns=[mean_sd_title, ALL_VISIT_DATES],
                     groupby=[stratify_analysis.STRATIFY_GROUP], overall=False, nonnormal=ALL_VISIT_DATES,
                     min_max=ALL_VISIT_DATES)
    return table


if __name__ == "__main__":
    main()
