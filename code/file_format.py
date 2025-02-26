#!/usr/bin/python3
"""
Hannah Takasuka, PhD student @UCSF Capra Lab, 2021-2024

\brief
    This script formats filenames for testing purposes, adding a datestamp and version to the end
    of the file for formal tests or overriding files for informal tests.
"""

import datetime
import os
import pathlib
from typing import Dict, List, Tuple

import pandas
import re

DEBUG = True
FILENAME = "test"
DATA_1 = {"one": 1, "two": 2, "three": 3}
DATA_2 = {"three": 3, "four": 4, "five": 5}
DATA_REPEAT = 5


def main() -> None:
    """Formulates file names, writes "hello world!" to a text file appends data to csv files."""
    filename_txt = format_filename(FILENAME, ".txt", DEBUG)
    with open(filename_txt, "w") as file:
        file.write("hello world!")
    filename_csv = format_filename(FILENAME, ".csv", DEBUG)
    for _ in range(DATA_REPEAT):
        record_dict_to_csv(filename_csv, DATA_1)
        record_dict_to_csv(filename_csv, DATA_2)
    add_filename_tag(filename_csv, "_tag")


def find_current_python_folder() -> pathlib.WindowsPath:
    return pathlib.Path(__file__).parent.resolve()


def find_recent_datetime_file_to_df(base_file: str, folder: str, nrows=None,
                                    index_col=None, exclude=None) -> pandas.DataFrame:
    filename = find_recent_datetime_file(base_file, folder, exclude)
    print(filename)
    df = pandas.read_csv(filename, nrows=nrows, index_col=index_col)
    return df


def find_recent_datetime_file(base_file: str, folder: str, exclude=None) -> str:
    """
    Finds the most recent file in the folder based on a datestamp in the filename
    For example, f("hello.txt") = "2050-01-01_hello.txt"
    """
    file_type, base_filename = separate_base_name_and_filetype(base_file)
    all_files = os.listdir(folder)
    filtered_files = []
    for file in all_files:
        if file_type in file and base_filename in file:
            if exclude:
                if exclude not in file:
                    filtered_files.append(f"{folder}\{file}")
            else:
                filtered_files.append(f"{folder}\{file}")
    latest_file = max(filtered_files, key=os.path.getctime)
    return latest_file


def find_filetype(filename: str) -> str:
    _, file_type = separate_base_name_and_filetype(filename)
    return file_type


def separate_base_name_and_filetype(combined_name: str) -> Tuple[str, str]:
    """Separate base name of a file from its type"""
    separator_index = combined_name.rfind(".")
    return combined_name[:separator_index], combined_name[separator_index:]


def format_filename_with_end_string(filename: str, end_string: str) -> str:
    base_name, file_type = separate_base_name_and_filetype(filename)
    return f"{base_name}{end_string}{file_type}"


def format_filename_with_date(filename: str, debug: bool) -> str:
    """Adds a date/time OR deletes and overrides the debug file name

    Creates filenames based on date/version for formal tests or basic filenames if in debug mode;
    Deletes previous debug files

    Args:
        filename: beginning name of the filename
        debug: If true, a simple overrideable filename is created. If false, a version-specific
               non-overrideable filename is created.

    Returns:
        filename: filename with date/version or basic debug filename
    """

    if debug:
        delete_file(filename)
    else:
        base_name, filetype = separate_base_name_and_filetype(filename)
        directory, base_name_no_directory = split_basename_directory(base_name)
        filename_nd = add_filename_version(base_name_no_directory, filetype)
        if directory:
            filename = f"{directory}/{filename_nd}"
        else:
            filename = f"{filename_nd}"
    return filename


def format_filename(base_name: str, filetype: str, debug: bool) -> str:
    """Adds a date/time OR deletes and overrides the debug file name

    Creates filenames based on date/version for formal tests or basic filenames if in debug mode;
    Deletes previous debug files

    Args:
        base_name: beginning name of the filename
        filetype: type such as ".csv", ".txt"
        debug: If true, a simple overrideable filename is created. If false, a version-specific
               non-overrideable filename is created.

    Returns:
        filename: filename with date/version or basic debug filename
    """
    if debug:
        filename = f"{base_name}{filetype}"
        delete_file(filename)
    else:
        directory, base_name_no_directory = split_basename_directory(base_name)
        filename_nd = add_filename_version(base_name_no_directory, filetype)
        if directory:
            filename = f"{directory}/{filename_nd}"
        else:
            filename = f"{filename_nd}"
    return filename


def delete_file(filename: str) -> None:
    """Deletes file if it exists or if its name exists as a part of a current file

    Args:
        filename: file to delete
    """
    if os.path.exists(filename):
        os.remove(filename)
    else:
        filetype_index = filename.rfind(".")
        filename_no_type = filename[:filetype_index]
        current_directory = [f for f in os.listdir(".") if os.path.isfile(f)]
        for file in current_directory:
            if filename_no_type in file:
                os.remove(file)


def add_filename_version(base_name: str, filetype: str) -> str:
    """Add a date stamp and version at the end of the string

    For example, if "2021-08-27_1_<base_name>.csv" already exists in the current directory, then
    this function will return "2021-08-27_2_<base_name>.csv"

    Args:
        base_name: beginning name of the filename
        filetype: type such as ".csv", ".txt"

    Returns:
        filename: "<date>_<version>_<base_name>.csv"
    """
    date_string = datetime.date.today().strftime("%Y-%m-%d")
    version = 1
    while True:
        proposed_filename = f"{date_string}_{version}_{base_name}"
        file_exists = False
        for file in os.listdir("."):
            if proposed_filename in file:
                file_exists = True
                break
        if file_exists:
            version += 1
        else:
            filename = f"{proposed_filename}{filetype}"
            break
    return filename


def split_basename_directory(base_name: str, separator="/") -> Tuple[str, str]:
    """Split basename and directory combined in 1 string to 2 different strings
    :param base_name: combined proposed filename and directory
    :param separator: character that separates subfolders in the string
    :returns: split name
    """
    if separator in base_name:
        split_directory = base_name.split(separator)
        base_name_no_directory = split_directory[-1]
        split_directory.pop(-1)
        directory = separator.join(split_directory)
        return directory, base_name_no_directory
    else:
        return "", base_name


def record_dict_to_csv(filename: str, *data: Dict) -> None:
    """Record dictionary data to a csv file

    Args:
        filename: filename to save data to
        data: data to save to the csv
    """

    concatenated_data = {}
    for element in data:
        concatenated_data.update(element)
    try:
        df_new = pandas.DataFrame(concatenated_data)
    except ValueError:
        df_new = pandas.DataFrame(concatenated_data, index=[0])
    if does_file_exist_with_data(filename):
        previous_data_df = pandas.read_csv(filename)
        df_new = pandas.concat([df_new, previous_data_df])
    df_new.to_csv(filename, index=False)


def does_file_exist_with_data(filename: str) -> bool:
    """Determines whether the file name exists and is not empty

    Args:
        filename: filename with applicable directory to check

    Returns:
        file_exists: whether the filename exists and is not empty
    """
    if os.path.exists(filename):
        return os.stat(filename).st_size > 0
    return False


def add_filename_tag(old_filename: str, new_tag: str) -> None:
    """Rename files with an added tag

    Args:
        old_filename: filename to be replaced
        new_tag: tag to add to old filename
    """
    filetype_index = old_filename.rfind(".")
    filename_no_type = old_filename[:filetype_index]
    file_type = old_filename[filetype_index:]
    new_filename = f"{filename_no_type}{new_tag}{file_type}"
    os.rename(old_filename, new_filename)


def read_csv_to_dict(filename: str) -> Dict:
    """Reads a csv to a dictionary.

    Args:
        filename: csv filename and directory, if applicable

    Returns:
        dictionary format of csv
    """
    df = pandas.read_csv(filename)
    data = {}
    for robot_object in df.columns:
        pose = []
        for joint_index in df.index:
            pose.append(df[robot_object][joint_index])
        try:
            robot_label = int(robot_object)
        except ValueError:
            robot_label = robot_object
        data[robot_label] = tuple(pose)
    return data


def round_dict(dictionary: Dict, round_digits=0) -> Dict:
    """Rounds the float values in a dictionary

    Args:
        dictionary: dictionary with float values
        round_digits: number of digits to round to
    """
    res = dict()
    for key, value in dictionary.items():
        if isinstance(value, float):
            if round_digits == 0:
                res[key] = round(value)
            else:
                res[key] = round(value, round_digits)
        else:
            res[key] = value
    return res


def keep_first_element_list(l: List) -> List:
    df = pandas.DataFrame(l, columns=['a'])
    return list(df.drop_duplicates(keep="first")['a'])


def find_first_number_in_string(s: str) -> int:
    all_numbers = re.findall(r'\d+', s)
    return int(all_numbers[0])


if __name__ == "__main__":
    main()
