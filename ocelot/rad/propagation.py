
def propagate(optics_line, dfl, optimize=False, dump=False):
    
    if optimize: #not implemented
        estimate_masks(optics_line, dfl)
        combine_elements(optics_line)
            
    for element in optics_line.sequence:
        element.mask.apply(dfl)
    return dfl

def estimate_masks(optics_line, dfl):#not implemented
    for element in optics_line.sequence:
        element.mask.get_mask(dfl)

def combine_elements(optics_line):#not implemented
    return optics_line