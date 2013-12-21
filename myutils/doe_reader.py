# This module is used to read csv files containing a specific set of experiments

import csv, os, sys

def index_in_unique_list(lst, a):
    ''' returns the index of where element a is present in lst 
        Assumption: lst is a list with ONLY unique elements
    '''
    if len(lst)!=len(set(lst)):
        print 'error: list is not unique'
        return 'error'
    else:
        tmp = [i for i, x in enumerate(lst) if x==a]
        res = tmp[0]
    return res

def read_experiment_field(doe_csv, header_name_of_exp_id, header_name):
    ''' this function takes a doe_csv fullname, opens the file,
    and returns a dictionary {'exp_id':'field_value'}
    Assumption: unique exp_id
    '''
    field_dict = {}
    with open(doe_csv,'rb') as csvfile:
        reader = csv.reader(csvfile)
        headers = reader.next()
        header_index_of_exp_id = index_in_unique_list(headers, header_name_of_exp_id)
        header_index_of_header_name = index_in_unique_list(headers, header_name)
        for row in reader:
            exp_id = row[header_index_of_exp_id]
            field_value = row[header_index_of_header_name]
            field_dict[exp_id] = field_value
    return field_dict

def create_experiment_fns(doe_csv, header_name_of_exp_id, in_dir, extension_name):
    ''' this function takes a doe_csv fullname, opens the file, 
    and returns a dictionary {'exp_id':'fullname_to_exp_file'}
    Assumption: unique exp_id, and all files in the same dir with the same extension
    '''
    field_dict = {}
    with open(doe_csv,'rb') as csvfile:
        reader = csv.reader(csvfile)
        headers = reader.next()
        header_index_of_exp_id = index_in_unique_list(headers, header_name_of_exp_id)
        for row in reader:
            exp_id = row[header_index_of_exp_id]
            field_dict[exp_id] = os.path.join(in_dir, exp_id+extension_name)
    return field_dict

def read_experiment_ids(doe_csv, header_name_of_exp_id):
    ''' this function takes a doe_csv fullname, opens the file,
    and creates a dictionary {'exp_id':'field_value'}
    Assumption: unique exp_id
    '''
    exp_ids = []
    with open(doe_csv,'rb') as csvfile:
        reader = csv.reader(csvfile)
        headers = reader.next()
        header_index_of_exp_id = index_in_unique_list(headers, header_name_of_exp_id)
        for row in reader:
            exp_id = row[header_index_of_exp_id]
            exp_ids.append(exp_id)
    return exp_ids
