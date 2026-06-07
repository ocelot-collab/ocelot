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


def check_value(value, value_ref, tolerance=1.0e-15, tolerance_type='relative', assert_info='', numerical_zero=1e-15):
    """Value with reference value check function

    if relative_tolerance='relative' then tolerance is relative (this is default value)
                                     if both numbers less than numerical zero then return None
    if relative_tolerance='absolute' then tolerance is absolute
    there is no tolerance for string values
    """

    if isinstance(value, str):
        if value == value_ref:
            return [None]
        else:
            return [assert_info + ' value is "' + value + '"\n reference value is "' + value_ref + '"\n\n']

    if tolerance_type == 'relative':
        if np.abs(value) < numerical_zero and np.abs(value_ref) < numerical_zero:
            return None

    if tolerance_type == 'relative':
        abs_value_ref = np.abs(value_ref)
    else:
        abs_value_ref = 1.0

    if np.abs(value - value_ref) <= tolerance * abs_value_ref:
        return None
    else:
        return assert_info + ' value is "' + str(value) + '"\n reference value is "' + str(value_ref) + '"\n tolerance is "' + str(tolerance) + '"\n tolerance type is "' + tolerance_type + '"\n\n'


def check_matrix(matrix, matrix_ref, tolerance=1.0e-15, tolerance_type='relative', assert_info='', numerical_zero=1e-15):
    """Matrix with reference matrix check function"""

    result = []
    for (index, x), (index_ref, x_ref) in zip(np.ndenumerate(matrix), np.ndenumerate(matrix_ref)):
        result.append(check_value(x, x_ref, tolerance, tolerance_type,
                                  assert_info=assert_info+' matrix element '+str(index)+'\n',
                                  numerical_zero=numerical_zero))
    
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
    """
    Legacy-compatible serializer:
    - 2D: returns a list of dicts { "0": v0, "1": v1, ... } for each ROW
    - 1D: returns a list of dicts with a single key "0": { "0": vi }  (so JSON stays a list)
    (Matches the old structure so existing refs keep working.)
    """
    a = np.asarray(np_array)

    # Enforce numeric floats for stability
    a = a.astype(float, copy=False)

    json_array = []

    if a.ndim == 1:
        # Represent each element as a single-entry dict so the outer structure stays a list
        for x in a:
            json_array.append({"0": float(x)})
        return json_array

    if a.ndim == 2:
        for row in a:
            line = {str(i): float(v) for i, v in enumerate(row)}
            json_array.append(line)
        return json_array

    # Optional: if higher-D appears, fall back to a meta format
    # (won't be used by your current tests, but keeps future door open)
    return {
        "__meta__": {"shape": a.shape, "dtype": str(a.dtype)},
        "data": a.ravel(order="C").tolist(),
    }


def json2numpy(json_array):
    """
    Legacy-compatible deserializer with ragged-row handling.
    - If dict with "__meta__": universal format -> reconstruct by shape.
    - If legacy list-of-dicts:
        * 1D case: each dict has a single key "0" -> returns shape (N,)
        * 2D case: rows may be ragged; we take the union of keys and fill missing with NaN.
    - Also tolerates rows given as plain lists (coerced to dicts with str(i) keys).
    Returns float dtype.
    """
    # Universal/meta branch (if you ever write such refs in future)
    if isinstance(json_array, dict) and "__meta__" in json_array:
        meta = json_array["__meta__"]
        shape = tuple(meta["shape"])
        data = np.array(json_array["data"], dtype=float)
        return data.reshape(shape, order="C")

    # Non-list: try direct numeric coercion
    if not isinstance(json_array, list):
        return np.array(json_array, dtype=float)

    # Empty list
    if len(json_array) == 0:
        return np.array([], dtype=float)

    # Normalize each row to a dict with string keys "0","1",...
    norm_rows = []
    line_lengths = []
    only_single_key = True

    for row in json_array:
        if isinstance(row, dict):
            # accept "0"/0 keys; standardize to string keys
            d = {}
            for k, v in row.items():
                kk = str(k)  # normalize key type
                d[kk] = v
            norm_rows.append(d)
            line_lengths.append(len(d))
            if len(d) != 1:
                only_single_key = False
        elif isinstance(row, (list, tuple, np.ndarray)):
            d = {str(i): row[i] for i in range(len(row))}
            norm_rows.append(d)
            line_lengths.append(len(d))
            if len(d) != 1:
                only_single_key = False
        else:
            # Scalar or unknown: treat as single value under key "0"
            norm_rows.append({"0": row})
            line_lengths.append(1)
            # only_single_key unchanged (len==1)

    # 1D legacy pattern: every row has exactly one key ("0")
    if all(L == 1 for L in line_lengths) and only_single_key:
        vals = []
        for d in norm_rows:
            # robustly read single key "0"
            vals.append(d.get("0", next(iter(d.values()))))
        return np.array(vals, dtype=float)

    # 2D case (possibly ragged): take union of keys
    # Convert keys to ints for stable numeric sort; ignore non-int keys safely
    all_keys = set()
    for d in norm_rows:
        for k in d.keys():
            try:
                all_keys.add(int(k))
            except Exception:
                # non-integer key: skip or handle separately; we skip to stay consistent with legacy numeric indexing
                pass

    if not all_keys:
        # No numeric keys found → try direct coercion
        return np.array(json_array, dtype=float)

    sorted_keys = sorted(all_keys)
    rows = []
    for d in norm_rows:
        row_vals = []
        for k in sorted_keys:
            v = d.get(str(k), np.nan)  # fill missing with NaN to keep rectangular
            row_vals.append(v)
        rows.append(row_vals)

    return np.array(rows, dtype=float)

def numpyBRT2json(np_array):
    return dict(enumerate(np_array.flatten()))


def json2numpyBRT(json_array):
    if len(json_array) <= 36:
        array = np.zeros((len(json_array)//6, 6))
        for i in json_array:
            array[int(i)//6, int(i)%6] = json_array[i]
    else:
        array = np.zeros((6, 6, 6))
        for i in json_array:
            x = json_array[i]
            i = int(i)
            array[i//36, (i - i//36 * 36)//6, (i - i//36 * 36)%6] = x
    return np.array(array, dtype=object)


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


def cmp_lists_of_mats(
    label: str,
    mats_cur,
    mats_ref_json,
    energies_list,
    results,
    tolerance: float = 1e-10,
    tolerance_type: str = "absolute",
    numerical_zero: float = 1e-15,
):
    """
    Compare lists of matrices (e.g., response matrices at several energies)
    using the existing `check_matrix` function and append messages to `results`.

    Parameters
    ----------
    label : str
        Short label (e.g. "Rx", "Ry", "Px", "Py") used in messages.
    mats_cur : list[np.ndarray]
        List of current matrices to test.
    mats_ref_json : list[dict]
        List of reference matrices serialized with numpy2json.
    energies_list : list[float]
        Energies (used for clear log messages).
    results : list
        Mutable list to which check messages (or None) are appended.
    tolerance : float, optional
        Numerical tolerance passed to `check_matrix`.
    tolerance_type : str, optional
        "absolute" or "relative" (as in `check_matrix`).
    numerical_zero : float, optional
        Threshold for treating values as zero.

    Notes
    -----
    * This helper is compatible with your `check_result()` framework
      — it appends `None` for pass and strings for fail.
    * It performs shape checks before calling `check_matrix`.
    """
    mats_ref = [json2numpy(j) for j in mats_ref_json]
    if len(mats_cur) != len(mats_ref):
        results.append(
            f"{label}: number of energy sets mismatch {len(mats_cur)} vs {len(mats_ref)}\n"
        )
        return

    for i, E in enumerate(energies_list):
        A = np.asarray(mats_cur[i])
        B = np.asarray(mats_ref[i])
        if A.shape != B.shape:
            results.append(
                f"{label} at E={E:g} GeV: shape mismatch {A.shape} vs {B.shape}\n"
            )
            continue

        msgs = check_matrix(
            A,
            B,
            tolerance=tolerance,
            tolerance_type=tolerance_type,
            assert_info=f"{label} at E={E:g} GeV - ",
            numerical_zero=numerical_zero,
        )
        results.extend(msgs)


def cmp_array(label: str, arr_cur, arr_ref_json, results,
              tolerance=1e-10, tolerance_type='absolute', numerical_zero=1e-15):
    """Compare one numeric array using your check_matrix + shape check."""
    A = np.asarray(arr_cur)
    B = json2numpy(arr_ref_json)

    if A.shape != B.shape:
        results.append(f"{label}: shape mismatch {A.shape} vs {B.shape}\n")
        return

    # check_matrix returns a list of None-or-message per element
    msgs = check_matrix(A, B, tolerance=tolerance, tolerance_type=tolerance_type,
                        assert_info=f"{label} - ", numerical_zero=numerical_zero)
    results.extend(msgs)