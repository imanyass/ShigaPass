import pandas as pd

def dict_to_str(kv_dict):
    result = ""
    for key, value in kv_dict.items():
        result += f"{key};{value}\n"
    return result

def write_out(content, outfpath, method="w"):
    with open(outfpath, method) as fout:
        fout.write(content)

def read_file(input_filepath):
    with open(input_filepath, "r") as fin:
        content = fin.read()
    return content

def get_val_from_file(csv_file, sort_column, value_column, delimiter=";", ascending=True, num_of_rows=1):
    """Get the value from the specified column in the top row that has been sorted via the sort column"""
    try:
        df = pd.read_csv(csv_file, delimiter=delimiter, header=None)
        sorted_df = df.sort_values(by=sort_column, ascending=ascending)
        row_values = list(sorted_df.iloc[:num_of_rows, value_column])
    except pd.errors.EmptyDataError:
        row_values = ""
    return row_values

def check_ascending(input_filepath, column_idx1, column_idx2):
    with open(input_filepath, "r") as fin:
        split_first_line = fin.readline().rstrip().split("\t")
        column1, column2 = split_first_line[column_idx1: column_idx2+1]
        if column1 > column2:
            return False
        return True
