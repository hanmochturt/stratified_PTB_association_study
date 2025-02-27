#!/usr/bin/python3
"""
Hannah Takasuka, rotation in the Capra Lab
January 2023
"""

import analysis
import file_format

import numpy as np
import pandas as pd

import stratified_plots
import sys
from typing import Dict, List, Set, Union

DEBUG = False
DEMOGRAPHICS_SAVE = True

FILENAME = f"{analysis.ADJUSTED_CRUDE_FILENAME}_phenotypes_pvalues.csv"
NOT_LABEL = "not_"
LESS_THAN_LABEL = f"less_than_"
THRESHOLD_AND_MORE_LABEL = "_and_up"
EARLY_PRETERM_THRESHOLD_WEEKS = 32
CONDITION_LABEL = f"below {EARLY_PRETERM_THRESHOLD_WEEKS}"
IND = "indicated"

# CSV labels
STRATIFY_COLUMN_LABELS = (analysis.GEST_WEEKS, analysis.PRETERM_TYPE)
STRATIFY_COLUMN_LABEL = STRATIFY_COLUMN_LABELS[1]  # TODO change this to not be hardcoded in functions
STRATIFY_DATA_VALUES = ["spontaneous", "PPROM", "PTL with TOCO and TERM"]
STRATIFY_GROUP = "stratify group"


def main() -> None:
    """Description"""
    if DEMOGRAPHICS_SAVE:
        stdout_origin = sys.stdout
        printout_filename_date = file_format.format_filename_with_date(
            f"early_vs_late_preterm_demographics.txt", DEBUG)
        sys.stdout = open(printout_filename_date, "w")

    cohort, birth_column_names, phecodes_key = analysis.obtain_data(DEBUG, analysis.DIRECTORY)
    cohort = analysis.define_preterm(cohort)

    if STRATIFY_COLUMN_LABEL == analysis.GEST_WEEKS:
        stratified_cohorts = stratify_early_late_preterm(cohort, add_term=True)
    elif STRATIFY_COLUMN_LABEL == analysis.PRETERM_TYPE:
        stratified_cohorts = stratify_spont_ind_preterm(cohort, add_term=True)

    analysis.print_demographics(cohort)
    print(STRATIFY_COLUMN_LABEL)
    print(analysis.obtain_distribution(STRATIFY_COLUMN_LABEL, cohort))
    # print(analysis.obtain_distribution(CONDITION_LABEL, cohort))

    for label, df in stratified_cohorts.items():
        preterm_dist = analysis.obtain_distribution(analysis.PRETERM, df)
        print(f"\n{label}: n="
              f"{df.shape[0] - (preterm_dist[False])}")
        if DEBUG:
            filename_date_version = file_format. \
                format_filename_with_date(f"{label}_stratified_data.csv", DEBUG)
            df.to_csv(filename_date_version, index=True)
        else:
            fit = analysis.measure_fit(df)
            phenotype_reg = analysis.convert_phecode_to_phenotype(fit, phecodes_key)
            phenotype_reg = phenotype_reg.sort_values(by=[analysis.P_VALUES])
            phenotype_reg.index.name = 'phenotype'
            filename_date_version = file_format.format_filename_with_date(f"{label}{FILENAME}",
                                                                          DEBUG)
            phenotype_reg.to_csv(filename_date_version, index=True)

    if DEMOGRAPHICS_SAVE:
        sys.stdout.close()
        sys.stdout = stdout_origin


def stratify_early_late_preterm(cohort: pd.DataFrame, add_term=False, term_separate=False) -> Dict:
    stratified_cohorts = stratify_cohort_number_threshold(cohort, analysis.GEST_WEEKS,
                                                          EARLY_PRETERM_THRESHOLD_WEEKS,
                                                          add_all_term=add_term)
    if term_separate:
        stratified_cohorts['Term'] = analysis.obtain_term_df(cohort)
    return stratified_cohorts


def stratify_spont_ind_preterm(cohort: pd.DataFrame, add_term=False, term_separate=False) -> Dict:
    cohort = drop_if_multiple_conditions(cohort, {analysis.PRETERM_TYPE: '0:No',
                                                  analysis.PRETERM: True})
    stratified_cohorts = {}
    if term_separate:
        stratified_cohorts['term'] = analysis.obtain_term_df(cohort)
        cohort = cohort.loc[cohort[analysis.PRETERM] == True]
    if term_separate or not add_term:
        cohort = analysis.obtain_preterm_df(cohort)
    preterm_stratified_cohorts = stratify_cohort_if_col_value(cohort, analysis.PRETERM_TYPE,
                                                              STRATIFY_DATA_VALUES,
                                                              add_all_term=add_term,
                                                              remove_full_term_from_in_group=True)
    stratified_cohorts.update(preterm_stratified_cohorts)
    stratified_cohorts[IND] = stratified_cohorts.pop("not_spontaneous")
    return stratified_cohorts


def add_spont_ind_term_column(cohort: pd.DataFrame) -> pd.DataFrame:
    stratified_cohort = stratify_spont_ind_preterm(cohort, term_separate=True)
    return move_dict_group_key_to_column(stratified_cohort)


def stratify_cohort_number_threshold(cohort: pd.DataFrame, column_title: str,
                                     column_threshold_value: Union[float, int],
                                     add_all_term=False,
                                     not_value_data=True) -> Union[Dict, pd.DataFrame]:
    in_group_indices = find_df_indices_if_element_number_threshold(cohort, column_title,
                                                                   column_threshold_value)
    out_group_indices = set(range(cohort.shape[0])) - in_group_indices
    if add_all_term:
        pre_term_group_indices = analysis.convert_true_f_list_to_true_indices_set(cohort[
                                                                                      analysis.PRETERM])
        term_group_indices = set(range(cohort.index.size, )) - pre_term_group_indices
        in_group_indices.update(term_group_indices)

    out_group_row_titles = cohort.index.values[list(out_group_indices)]
    in_group_row_titles = set(cohort.index.values[list(in_group_indices)])
    cohort_in = cohort.filter(in_group_row_titles, axis='index')

    cohort_null = cohort.filter(out_group_row_titles, axis='index')
    if not_value_data:
        return {f"less_than_{column_threshold_value}": cohort_in,
                f"{column_threshold_value}_and_up": cohort_null}
    else:
        return cohort_in


def stratify_cohort_if_col_value(cohort: pd.DataFrame, column_title: str, column_values: List,
                                 remove_full_term_from_in_group=False, add_all_term=False,
                                 not_value_data=True) -> \
        Union[
            Dict, pd.DataFrame]:
    not_applicable_indices = find_df_indices_if_element_in_column(cohort, column_title, ["NA"])
    in_group_indices = find_df_indices_if_element_in_column(cohort, column_title, column_values) - \
                       not_applicable_indices
    # TODO not applicable should be for 37+ weeks but not indicated/sponteanous
    out_group_indices = set(range(cohort.shape[0])) - in_group_indices - not_applicable_indices

    if add_all_term or remove_full_term_from_in_group:
        pre_term_group_indices = analysis.convert_true_f_list_to_true_indices_set(cohort[
                                                                                      analysis.PRETERM])
        term_group_indices = set(range(cohort.index.size, )) - pre_term_group_indices
        row_titles = cohort.index.values[list(not_applicable_indices & pre_term_group_indices)]
        print(f"{len(not_applicable_indices & pre_term_group_indices)} preterm & not applicable "
              f"{column_title} dropped from "
              f"both stratify groups. Their IDs: {row_titles}\n")
    if remove_full_term_from_in_group:
        term_and_spontaneous = in_group_indices.intersection(term_group_indices)
        term_and_spontaneous_row_titles = cohort.index.values[list(term_and_spontaneous)]
        print(f"There are {len(term_and_spontaneous)} term and {column_values[0]} patients. Their"
              f" IDs: {term_and_spontaneous_row_titles}\n")
    if add_all_term:
        in_group_indices.update(term_group_indices)

    out_group_row_titles = cohort.index.values[list(out_group_indices)]
    in_group_row_titles = cohort.index.values[list(in_group_indices)]
    cohort_in = cohort.filter(in_group_row_titles, axis='index')
    cohort_null = cohort.filter(out_group_row_titles, axis='index')
    if not_value_data:
        return {column_values[0]: cohort_in, f"{NOT_LABEL}{column_values[0]}": cohort_null}
    else:
        return cohort_in


def find_df_indices_if_element_number_threshold(df: pd.DataFrame, column_title: str,
                                                column_threshold_value: Union[float, int]) -> Set:
    in_group_true_false = df[column_title] < column_threshold_value
    in_group_true_false = in_group_true_false.to_numpy()
    in_group_indices = analysis.convert_true_f_list_to_true_indices_set(list(in_group_true_false))
    return in_group_indices


def find_df_indices_if_element_in_column(df: pd.DataFrame, column_title: str,
                                         column_values: List) -> Set:
    # TODO: Clean up this code
    all_in_group_true_false = np.full((df.index.size,), False)
    for column_value in column_values:
        if isinstance(column_value, str):
            try:
                in_group_true_false = df[column_title].str.contains(column_value).to_numpy()
            except AttributeError:
                if column_value == True:
                    in_group_true_false = df[column_title].to_numpy()
                else:
                    if column_value == False:
                        raise ValueError("This code probably doesn't work for column_value = False."
                                         "Feel free to delete this raised error and see if it works correctly.")
                    in_group_true_false = (df[column_title] == column_value).to_numpy()
        else:
            if column_value == True:
                in_group_true_false = df[column_title].to_numpy()
            else:
                if column_value == False:
                    raise ValueError("This code probably doesn't work for column_value = False."
                                     "Feel free to delete this raised error and see if it works correctly.")
                in_group_true_false = (df[column_title] == column_value).to_numpy()
        all_in_group_true_false = np.logical_or(all_in_group_true_false, in_group_true_false)
    in_group_indices = analysis.convert_true_f_list_to_true_indices_set(list(
        all_in_group_true_false))
    return in_group_indices


def drop_if_multiple_conditions(df: pd.DataFrame, condition: Dict) -> pd.DataFrame:
    """drop rows from dataframe if they have particular values in their columns

    :param condition: column title (key), column value (value) to be removed"""
    in_group_row_titles = find_indices_if_multiple_conditions(df, condition)
    return df.filter(in_group_row_titles, axis='index')


def find_indices_if_multiple_conditions(df: pd.DataFrame, condition: Dict) -> List:
    condition_indcies = None
    for column_title, column_value in condition.items():
        if not isinstance(column_value, (list, tuple)):
            column_value = [column_value]
        col_values_indices = find_df_indices_if_element_in_column(df, column_title, column_value)
        if not condition_indcies:
            condition_indcies = col_values_indices
        else:
            condition_indcies = col_values_indices & condition_indcies
    print(
        f"\n{len(condition_indcies)} rows that meet the conditions {condition} were dropped. Their IDs: {df.index.values[list(condition_indcies)]}")
    in_group_indices = set(range(df.shape[0])) - condition_indcies
    in_group_row_titles = df.index.values[list(in_group_indices)]
    return in_group_row_titles


def move_dict_group_key_to_column(cohorts: Dict) -> pd.DataFrame:
    dfs_with_group_titles_column = []
    for group_title, group_cohort in cohorts.items():
        group_cohort[STRATIFY_GROUP] = group_title
        dfs_with_group_titles_column.append(group_cohort)
    return pd.concat(dfs_with_group_titles_column)


def move_dict_group_key_to_column(cohorts: Dict) -> pd.DataFrame:
    dfs_with_group_titles_column = []
    for group_title, group_cohort in cohorts.items():
        print(group_title, len(group_cohort))
        group_cohort[STRATIFY_GROUP] = group_title
        dfs_with_group_titles_column.append(group_cohort)
    return pd.concat(dfs_with_group_titles_column)


if __name__ == "__main__":
    main()
