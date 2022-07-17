import csv


def read_csv(file_name):
    """This function reads the mutation data in csv format and returns a list."""
    mut_list = []
    with open(file_name, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mut_list.append(row)
    return mut_list