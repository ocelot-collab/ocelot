"""Common functions and parameters"""

import numpy as np
import json
import pytest 

"""Tolerance parameter
 in the case of using relative tolerance
 absolute(value - reference_value) <= TOL * absolute(reference_value)

 in the case of using absotute tolerance
 absolute(value - reference_value) <= TOL
"""
TOL = 1.0e-7


def add_info(data):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(data + '\n')
    f.close()


def check_result(data):
    
    result = True
    for r in data:
        if r is None:
            if result:
                info = ' [PASSED]\n'

        else:
            if result:
                info = ' [FAILED]\n'

            info += r
            result *= False

    add_info(info)
    return result


def check_value(value, value_ref, tolerance=1.0e-15, tolerance_type='relative', assert_info=''):
    """Value with reference value check function

    if relative_tolerance='relative' then tolerance is relative (this is default value)
    if relative_tolerance='absotute' then tolerance is absotute
    there is no tolerance for string values
    """

    if isinstance(value, str):
        if value == value_ref:
            return [None]
        else:
            return [assert_info + ' value is "' + value + '"\n reference value is "' + value_ref + '"\n\n']

    if tolerance_type == 'relative':
        abs_value_ref = np.abs(value_ref)
    else:
        abs_value_ref = 1.0

    if np.abs(value - value_ref) <= tolerance * abs_value_ref:
        return None
    else:
        return assert_info + ' value is "' + str(value) + '"\n reference value is "' + str(value_ref) + '"\n tolerance is "' + str(tolerance) + '"\n tolerance type is "' + tolerance_type + '"\n\n'


def check_matrix(matrix, matrix_ref, tolerance=1.0e-15, tolerance_type='relative', assert_info=''):
    """Matrix with reference matrix check function"""

    delta_matrix = np.abs(matrix - matrix_ref)
    if tolerance_type == 'relative':
        tol_matrix = tolerance * np.abs(matrix_ref)
    else:
        tol_matrix = tolerance

    result = []
    for (index, x), (index_ref, x_ref) in zip(np.ndenumerate(matrix), np.ndenumerate(matrix_ref)):
        result.append(check_value(x, x_ref, tolerance, tolerance_type, assert_info=assert_info+' matrix element '+str(index)+'\n'))
    
    return result


def check_dict(dict_t, dict_ref, tolerance=1.0e-15, tolerance_type='relative', assert_info=''):
    """Dictionary with reference dictionary check function"""
    
    if len(dict_t) != len(dict_ref):
        return [assert_info + '\n dict len is "' + str(len(dict_t)) + '"\n reference dict len is "' + str(len(dict_ref)) + '"\n']

    result = []
    line = -1
    for elem, elem_ref in zip(dict_t, dict_ref):
        line += 1
        common_keys = set(elem.keys()) & set(elem_ref.keys())
        for key in common_keys:

            if key == 'id':
                continue
            
            if isinstance(elem[key], list):
                result += check_matrix(np.asarray(elem[key]), np.asarray(elem_ref[key]),
                                       tolerance, tolerance_type, assert_info=assert_info+
                                                                              ' line is '+str(line)+ "/" + str(len(dict_t)) +' key is '+str(key)+'\n')
                continue

            if isinstance(elem[key], dict):
                result += check_dict([elem[key]], [elem_ref[key]],
                                     tolerance, tolerance_type, assert_info=assert_info+
                                                                            ' line is '+str(line)+ "/" + str(len(dict_t)) + ' key is '+str(key)+'\n')
                continue

            result.append(check_value(elem[key], elem_ref[key],
                                      tolerance, tolerance_type, assert_info=assert_info+
                                                                             'line is '+str(line)+ "/" + str(len(dict_t)) +' key is '+str(key)+'\n'))
    
    return result


def numpy2json(np_array):
    
    json_array = []

    for column in np_array:
        line = dict(enumerate(column))
        json_array.append(line)
    
    return json_array


def json2numpy(json_array):

    array = []
    
    for line in json_array:
        tmp = []
        for i in range(len(line)):
            tmp.append(line[str(i)])
        array.append(tmp)

    return np.array(array)


def obj2dict(object, unpack=None):
    
    json_obj = []

    for elem in object:
        tmp = {}
        for key, value in elem.__dict__.items():

            if isinstance(value, (float, int)):
                tmp[key] = value
                continue
            
            if isinstance(value, list):
                tmp[key] = value
                continue

            if unpack is not None and key in unpack:
                tmp_p = {}
                for key_p, value_p in value.__dict__.items():
                    tmp_p[key_p] = value_p
                tmp[key] = tmp_p
                continue
                
        json_obj.append(tmp)

    return json_obj


def json_read(filename):
    return json.load(open(filename))
    
    
def json_save(data, filename):
    with open(filename, "w") as text_file:
        text_file.write(json.dumps(data))


