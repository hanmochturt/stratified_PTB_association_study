#!/usr/bin/python3
"""
Hannah Takasuka, rotation in the Capra Lab
January 2023

parts in Jean Costello's R code that I didn't incorporate:
* cohorts split by all deliveries, first pregnancies, and those with pregnancy complications
diagnosis codes
* proportion difference because this doesn't seem to be reflected in the manhattan plot
"""
import decimal
import math
from math import log10, floor
import pathlib
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as smm
import sys
from typing import Dict, List, Set, Tuple, Union

from itertools import compress
import numpy as np
import pandas as pd
from patsy import dmatrix
from tqdm import tqdm

import file_format

DEBUG = True  # only use first 100 entries to make the code run faster
ADJUSTED_ODDS = True
DEMOGRAPHICS_SAVE = False

ADJUSTED_CRUDE_FILENAME = {True: "_adjusted", False: "_crude"}[ADJUSTED_ODDS]
FILENAME = f"not_stratified_phenotypes_pvalues{ADJUSTED_CRUDE_FILENAME}.csv"

# dataframe labels (can be edited)
PRETERM = "preterm"
UNKNOWN = "unknown"
EDU_SIMP = "MATEDUC_simplified"

# output CSV labels (can be edited)
P_VALUES = "pvalues"
P_VALUES_FDR = "pvalues_fdr_corrected"
OR_COEFF = "or"
OR_CI_LOWER = "or_lower"
OR_CI_UPPER = "or_upper"
N_PRETERM = "n_preterm"
N_TOTAL = "n_total"
PHENOTYPE_CATEGORY = "phe_category"

# input CSV labels
GEST_WEEKS = "GAWEEKS"
PRETERM_TYPE = "PRETERM"
RACE = "race9c"
INS_PRIVATE = "INSPRIVATE"
PHECODE = "phecode"
PHENOTYPE = "phenotype"
PHE_CATEGORY = "category"
AGE = "MATAGE"
MAT_EDU = "MATEDUC"
ID_BIRTHS = 'person_id_mom'
ID_INF = 'person_id_inf'

# existing folder names
BASE_DIRECTORY_KEY = 'base'
DEF_FOLDER = 'defs'
DATA_FOLDER = 'data'
INTERMEDIATE_FOLDER = 'intermediate'
IMAGE_FOLDER = 'images'
OUTPUT_FOLDER = 'output'

# existing file names
CONDITIONS = 'conditions_binary'
BIRTHS = 'births'
PHECODES = 'phecode_definitions1.2.csv'
DIRECTORY = r"//ars-data-01.sde.net.ucsf.edu/MyResearchShared/sirotam1_shared/Birth_IRB1722929/projects/ptb_diagnoses"


def main() -> None:
    print(DIRECTORY, "base dir begin")
    """Description"""
    if DEMOGRAPHICS_SAVE:
        stdout_origin = sys.stdout
        printout_filename_date = file_format.format_filename_with_date(f"all_demographics.txt", DEBUG)
        sys.stdout = open(printout_filename_date, "w")
    cohort, birth_column_names, phecodes_key = obtain_data(DEBUG, DIRECTORY,
                                                           cond_file=f"{DIRECTORY}/data/2023-07-18_conditions_binary.csv",)# births_file=f"{DIRECTORY}/data/2023-04-27_births.csv")
    cohort = define_preterm(cohort)
    print_demographics(cohort)

    fit = measure_fit(cohort)

    phenotype_reg = convert_phecode_to_phenotype(fit, phecodes_key)
    phenotype_reg = phenotype_reg.sort_values(by=[P_VALUES])
    filename_date_version = file_format.format_filename_with_date(FILENAME, DEBUG)
    phenotype_reg.to_csv(filename_date_version, index=True)

    if DEMOGRAPHICS_SAVE:
        sys.stdout.close()
        sys.stdout = stdout_origin


def obtain_data(debug=False, base_directory=None, cond_file=None,
                births_file=None, reduce_phecodes=False) -> \
        pd.DataFrame:
    """ open relevant csv files to Pandas Dataframes based on directories defined in the constants
    :param debug: whether a shortened dataframe should be obtained to make later data analysis
    quicker
    :returns: dataframes of data extracted from the directories
    """
    directories = set_directories(base_directory)
    if debug:
        n_rows = 100
    else:
        n_rows = None
    if not cond_file:
        cond_df = file_format.find_recent_datetime_file_to_df(CONDITIONS, directories[
            DATA_FOLDER], n_rows)
    else:
        cond_df = pd.read_csv(cond_file, nrows=n_rows)
    births_df = obtain_birth_data(DIRECTORY, births_file, n_rows)
    phecodes_df = pd.read_csv(f"{directories[DEF_FOLDER]}/{PHECODES}")
    birth_column_names = births_df.columns.tolist()
    birth_column_names.remove(ID_BIRTHS)
    cohort = merge_births_conds(births_df, cond_df)
    if reduce_phecodes:
        cohort = reduce_phecodes_to_first_decimal(cohort)
    for column in "INSPRIVATE", "MATEDUC":
        cohort[column] = cohort[column].fillna("unknown")
    return cohort, birth_column_names, phecodes_df


def obtain_birth_data(base_directory=None, file=None, n_rows=None) -> pd.DataFrame:
    directories = set_directories(base_directory)
    if not file:
        births_df = file_format.find_recent_datetime_file_to_df(BIRTHS, directories[
            DATA_FOLDER])
    else:
        births_df = pd.read_csv(file, nrows=n_rows)
    births_df = define_preterm(births_df)
    return births_df


def set_directories(base_directory=None) -> Tuple:
    """create relevant directory names based on defaulting to the parent of the current folder or
    what is user-defined

    :param base_directory: string of base directory, or none if running the directory of the script
    :returns: directories of multiple folders
    """
    if base_directory:
        filepath = pathlib.Path(base_directory)
    else:
        filepath = pathlib.Path(__file__).parent.resolve()
        print(filepath)
    directories = {BASE_DIRECTORY_KEY: str(filepath)}
    for folder in (DEF_FOLDER, DATA_FOLDER, INTERMEDIATE_FOLDER, IMAGE_FOLDER, OUTPUT_FOLDER):
        directories[folder] = str(filepath / folder)
    return directories


def define_preterm(births: pd.DataFrame) -> pd.DataFrame:
    """Add a column to a dataframe that is true/false for whether gestational age is below 37
    """
    preterm_binary = births[GEST_WEEKS] < 37
    return pd.concat((preterm_binary.rename(PRETERM), births), axis=1)


def merge_births_conds(births: pd.DataFrame, conds: pd.DataFrame) -> pd.DataFrame:
    try:
        print("Used infant ID to link births and conditions.")
        conds = conds.set_index('person_id_inf')
        births = births.set_index('person_id_inf')
    except KeyError:
        print("Used maternal ID to link births and conditions.")
        conds = conds.set_index('person_id')
        births = births.set_index('person_id_mom')
    try:
        cohort = pd.concat([births, conds], axis=1, join='inner')
    except pd.errors.InvalidIndexError:
        conds = conds[~conds.index.duplicated(keep='first')]
        births = births[~births.index.duplicated(keep='first')]
        cohort = pd.concat([births, conds], axis=1, join='inner')
    return cohort


def convert_phecode_str_to_float(phecode: str) -> float:
    return float(str(phecode.replace("P", "")))


def truncate(number, digits=1) -> float:
    # Improve accuracy with floating point operations, to avoid truncate(16.4, 2) = 16.39 or truncate(-1.13, 2) = -1.12
    nbDecimals = len(str(number).split('.')[1])
    if nbDecimals <= digits:
        return number
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


def find_pairs_common_first_decimal(float_sequence: List):
    float_sequence.sort()
    same_first_decimal = []
    unique_first_decimal = []
    previous_paired = False
    for index, phecode in enumerate(float_sequence[:-1]):
        if truncate(float_sequence[index+1], 1) == truncate(phecode, 1):
            if previous_paired:
                same_first_decimal[-1].append(float_sequence[index+1])
            else:
                same_first_decimal.append([phecode, float_sequence[index+1]])
            previous_paired = True
        else:
            unique_first_decimal.append(phecode)
            previous_paired = False
    if truncate(float_sequence[-2]) == truncate(float_sequence[-1]):
        unique_first_decimal.append(float_sequence[-1])
    return same_first_decimal, unique_first_decimal


def convert_float_to_str_phecode(number: float) -> str:
    return f"P{number}"


def reduce_phecodes_to_first_decimal(cohort: pd.DataFrame) -> pd.DataFrame:
    phecode_cols = get_phecode_columns(cohort)
    phecode_float = list(map(convert_phecode_str_to_float, phecode_cols))
    str_to_float_key = dict(zip(phecode_cols, phecode_float))
    float_to_str_key = dict(zip(phecode_float, phecode_cols))
    cohort_float = cohort.rename(columns = str_to_float_key)
    same_first_decimal_pairs, unique_decimal_phecodes = find_pairs_common_first_decimal(phecode_float)
    phecode_truncated_float_cols = []
    for same_first_decimal_pair in same_first_decimal_pairs:
        cohort_float[truncate(same_first_decimal_pair[0])] = cohort_float[same_first_decimal_pair].sum(axis=1).astype(bool)
        cohort_float[truncate(same_first_decimal_pair[0])] = cohort_float[truncate(same_first_decimal_pair[0])].astype(int)
        phecode_truncated_float_cols.append(truncate(same_first_decimal_pair[0]))
        if truncate(same_first_decimal_pair[0]) == same_first_decimal_pair[0]:
            same_first_decimal_pair.pop(0)
        cohort_float = cohort_float.drop(columns=same_first_decimal_pair)
    phecode_truncated_str_cols = cohort_float[phecode_truncated_float_cols].rename(columns=convert_float_to_str_phecode)
    unique_decimal_phecodes_dict = dict((key, float_to_str_key[key]) for key in unique_decimal_phecodes)
    string_unique_decimal_phecodes = unique_decimal_phecodes_dict.values()
    phecode_unique_decimal_cols = cohort[string_unique_decimal_phecodes]
    cohort_no_phecodes = cohort.drop(columns=phecode_cols)
    cohort_str_tuncated = pd.concat([cohort_no_phecodes, phecode_truncated_str_cols, phecode_unique_decimal_cols], axis=1)
    return cohort_str_tuncated[~cohort_str_tuncated.index.duplicated(keep='first')]


def obtain_distribution(df_column: str, cohort: pd.DataFrame) -> Dict:
    series = cohort[df_column]
    series_na_label = series.fillna("NaN")
    dist = series_na_label.value_counts()
    dist_keys = list(dist.keys())
    try:
        dist_keys.sort()
        sorted_dict = {i: dist[i] for i in dist_keys}
        return sorted_dict
    except TypeError:
        return dist


def print_demographics(cohort: pd.DataFrame) -> None:
    cohort = add_maternal_edu_simplified_column(cohort)
    for column in (RACE, PRETERM, INS_PRIVATE, EDU_SIMP):
        print(column)
        print(obtain_distribution(column, cohort), '\n')


def measure_fit(cohort: pd.DataFrame) -> pd.DataFrame:
    if ADJUSTED_ODDS:
        cohort = insert_spline_in_df(cohort, AGE)
        cohort = add_maternal_edu_simplified_column(cohort) #TODO: check if this adjustment works
    formula_inputs_sum = create_fit_formula_sum_inputs(cohort, ADJUSTED_ODDS)
    or_data = pd.DataFrame()
    phecodes = get_phecode_columns(cohort)
    for phecode in tqdm(phecodes):
        if len(list(obtain_distribution(phecode, cohort).keys())) > 1:
            fit_data = measure_fit_phecode(cohort, phecode, formula_inputs_sum)
            or_data = pd.concat([or_data, fit_data])
    or_data[P_VALUES_FDR] = calculate_p_value_fdr_correction(list(or_data[P_VALUES]))
    return or_data


def insert_spline_in_df(cohort: pd.DataFrame, column_title: str, df=3) -> pd.DataFrame:
    spline = create_spline(cohort, column_title, df)
    data_with_spline = pd.concat([spline, cohort], axis='columns')
    return data_with_spline.drop(columns=column_title)


def create_spline(cohort: pd.DataFrame, column_title: str, df=3) -> pd.DataFrame:
    print((f"{cohort[column_title].isnull().values.sum()} NA value from the {column_title} "
           f"column removed."))
    cohort = cohort[cohort[column_title].notna()]
    spline_matrix = dmatrix(f"cr(x, df={df}) - 1", {"x": cohort[column_title]},
                            return_type='dataframe')
    for i, column in enumerate(spline_matrix.columns):
        spline_matrix = spline_matrix.rename(columns={column: f"{column_title}_{i}"})
    return spline_matrix


def measure_fit_phecode(cohort, phecode, formula_inputs_sum=None, adjusted_odds=ADJUSTED_ODDS):
    original_phecode = phecode
    if "." in str(original_phecode):
        phecode_nickname = original_phecode.replace(".", "decimal")
        cohort = cohort.rename(columns={original_phecode: phecode_nickname})
        phecode = phecode_nickname
    if not formula_inputs_sum:
        formula_inputs_sum = create_fit_formula_sum_inputs(cohort, adjusted_odds)
    formula = f"{PRETERM} ~ {formula_inputs_sum} + {phecode}"
    model = smf.glm(formula=formula,
                    data=cohort,
                    family=sm.families.Binomial())
    result = model.fit()
    coeff_upper_lower = [result.params[phecode],
                         result.conf_int(alpha=0.05, cols=None)[0][phecode],
                         result.conf_int(alpha=0.05, cols=None)[1][phecode]]
    or_upper_lower = []
    for value in coeff_upper_lower:
        try:
            odds_ratio = math.exp(value)
        except OverflowError:
            odds_ratio = 1E10
        if odds_ratio == 0:
            or_upper_lower.append(1E-10)  # avoid zero division error
        else:
            or_upper_lower.append(odds_ratio)
    # TODO: investigate why odds ratios are the inverse of what they should be
    try:
        index_name = [float(str(original_phecode.replace("P", "")))]
    except ValueError:
        index_name = [original_phecode]
    fit_data = pd.DataFrame({OR_COEFF: 1 / or_upper_lower[0],
                             OR_CI_LOWER: 1 / or_upper_lower[2],
                             OR_CI_UPPER: 1 / or_upper_lower[1],
                             P_VALUES: result.pvalues[phecode],
                             N_PRETERM: np.sum(obtain_preterm_df(cohort)[
                                                   phecode].to_numpy()),
                             N_TOTAL: np.sum(cohort[phecode].to_numpy())},
                            index=index_name)
    return fit_data


def create_fit_formula_sum_inputs(cohort: pd.DataFrame, adjusted_odds: bool, age=True):
    if adjusted_odds:
        if age:
            formula_inputs = [i for i in list(cohort.columns) if AGE in i]
        else:
            formula_inputs = []
        formula_inputs.append(EDU_SIMP)
        formula_inputs.append(INS_PRIVATE)
        formula_inputs.append(RACE)
        formula_inputs_sum = " + ".join(formula_inputs)
    else:
        formula_inputs_sum = ""
    return formula_inputs_sum


def remove_na_from_dict(data: Dict) -> Dict:
    for key, values_list in data.items():
        if np.isnan(values_list).any():
            values_numpy = np.array(values_list)
            data[key] = list(values_numpy[~np.isnan(values_numpy)])
            print(f"{len([~np.isnan(values_numpy)])} NA value removed from {key}")
    return data


def obtain_preterm_df(cohort: pd.DataFrame) -> pd.DataFrame:
    pre_term_group_indices = convert_true_f_list_to_true_indices_set(cohort[PRETERM])
    return obtain_df_from_number_indices(cohort, pre_term_group_indices)


def obtain_term_df(cohort: pd.DataFrame) -> pd.DataFrame:
    pre_term_group_indices = convert_true_f_list_to_true_indices_set(cohort[PRETERM])
    term_group_indices = set(range(cohort.index.size, )) - pre_term_group_indices
    return obtain_df_from_number_indices(cohort, term_group_indices)


def obtain_df_from_number_indices(cohort: pd.DataFrame, indices: Union[Set, List, Tuple]) -> \
        pd.DataFrame:
    in_group_row_titles = cohort.index.values[list(indices)]
    return cohort.filter(in_group_row_titles, axis='index')


def convert_true_f_list_to_true_indices_set(true_false_list) -> Set:
    if not isinstance(true_false_list, list):
        true_false_list = list(true_false_list)
    return set(compress(range(len(true_false_list)), true_false_list))


def convert_phecode_to_phenotype(phecode_regression: pd.DataFrame, phecodes_key: pd.DataFrame) ->\
        pd.DataFrame:
    phecode_pheno_dict = pd.Series(phecodes_key[PHENOTYPE].values, index=phecodes_key[
        PHECODE])
    phecode_category_dict = pd.Series(phecodes_key[PHE_CATEGORY].values, index=phecodes_key[
        PHECODE]).to_dict()
    phecode_regression[PHENOTYPE_CATEGORY] = phecode_regression.index.map(
        phecode_category_dict)
    phenotype_regression = phecode_regression.rename(index=phecode_pheno_dict)
    return phenotype_regression


def convert_phenotype_to_phecode(phenotype, phecodes_key):
    phecodes_key_dict = dict(zip(phecodes_key[PHENOTYPE].values, list(phecodes_key[PHECODE])))
    return phecodes_key_dict[phenotype]


def calculate_p_value_fdr_correction(p_values: Union[List, Tuple]) -> Tuple:
    nan_indices = np.where(np.isnan(p_values))[0]
    fdr_pval = np.empty((len(p_values)))
    fdr_pval[nan_indices] = np.nan
    pval_indices = set(range(len(p_values))) - set(nan_indices)

    pval_no_nan = [x for x in p_values if x == x]
    fdr_pval_no_nan = smm.fdrcorrection(pval_no_nan)[1]
    for index, pval_index in enumerate(pval_indices):
        fdr_pval[pval_index] = fdr_pval_no_nan[index]
    return fdr_pval


def round_sig_figs(number, num_figs=1):
    if number == 0:
        new_number = number
    elif math.isnan(number):
        new_number = math.nan
    else:
        new_number = round(number, num_figs-1 - int(floor(log10(abs(number)))))
    return new_number


def round_sig_figs_pandas(df, num_figs=1):
    if isinstance(df, pd.DataFrame):
        df = df.applymap(lambda x: round_sig_figs(x, num_figs))
        df = df.applymap(lambda x: format_number(x))
    elif isinstance(df, pd.Series):
        df = df.map(lambda x: round_sig_figs(x, num_figs))
        df = df.map(lambda x: format_number(x))
    else:
        print("type is not pandas dataframe or series. values were not rounded.")
    return df


def format_number(num):
    """remove trailing zeros from numbers after decimals"""
    if num == 0:
        return 0
    else:
        num = str(num)
    try:
        dec = decimal.Decimal(num)
    except:
        return num
    tup = dec.as_tuple()
    delta = len(tup.digits) + tup.exponent
    digits = ''.join(str(d) for d in tup.digits)
    if delta <= 0:
        zeros = abs(tup.exponent) - len(tup.digits)
        val = '0.' + ('0'*zeros) + digits
    else:
        val = digits[:delta] + ('0'*tup.exponent) + '.' + digits[delta:]
    val = val.rstrip('0')
    if val[-1] == '.':
        val = val[:-1]
    if tup.sign:
        return '-' + val
    return val


def get_phecode_columns(df: pd.DataFrame) -> List:
    phecode_labels = []
    for column_label in df.columns:
        if "P" in column_label[0] and column_label[1].isnumeric():
            phecode_labels.append(column_label)
    return phecode_labels


def add_maternal_edu_simplified_column(cohort: pd.DataFrame) -> pd.DataFrame:
    edu_map = make_grade_level_edu_map(cohort)
    cohort_edu = cohort[MAT_EDU].map(edu_map.get)
    cohort[EDU_SIMP] = pd.Categorical(cohort_edu, categories=set(edu_map.values()))
    return cohort


def obtain_int_from_str(string: str):
    num_string = ""
    for character in string:
        if character.isdigit():
            num_string += character
    if num_string:
        return int(num_string)
    else:
        print(f"{string} value could not be reduced to an integer")
        return string


def make_grade_level_edu_map(cohort) -> Dict:
    edu = obtain_distribution(MAT_EDU, cohort)
    edu_map = {}
    for string in edu.keys():
        grade_level_int = obtain_int_from_str(string)
        try:
            if grade_level_int < 12:
                edu_map[string] = "less_than_12th_grade"
            elif grade_level_int == 12:
                edu_map[string] = "12th_grade"
            elif grade_level_int > 12:
                edu_map[string] = "college"
        except TypeError:
            edu_map[string] = string
    return edu_map


if __name__ == "__main__":
    main()
