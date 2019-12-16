
def propagate(optics_line, dfl, optimize=True, dump=False):
    if optimize:
        estimate_masks(optics_line, dfl)
        combine_elements(optics_line)

    for element in optics_line.sequence:

        element.mask.apply(dfl)
    return dfl

def estimate_masks(optics_line, dfl):
    for element in optics_line.sequence:
        element.mask.get_mask(dfl)

def combine_elements(optics_line):
    return optics_line