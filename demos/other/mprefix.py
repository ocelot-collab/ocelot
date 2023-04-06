# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: sserkez
"""

import numpy as np
import time
# from ocelot.common.math_op import mprefix


def mprefix(value, order_bias=0.05, apply=1):
    '''
    estimate metric prefix for a scalar or np.array
    accepts floating point number
    order_bias is an ammount of order of magnitudes to round up to. Allows to avoid values like 995 instead of 0.995k. 1/3 corresponds to ...980,990,0.10,0.11...
    returns scaled values (float or numpy arrayof floats), its prefix (string), and applied order of magnitude (should it need to be reused)
        '''
    prefix_range = (-30, 30)
    prefixes = {-30:'q', -27:'r', -24:'y', -21:'z', -18:'a', -15:'f', -12:'p', -9:'n', -6:r'$\mu$', -3:'m', 0:'', 3:'k', 6:'M', 9:"G", 12:'T', 15:'P', 18:'E', 21:'Z', 24:'Y', 27:'R', 30:'Q'}
        
    if not apply:
        return (value, '', 0)
    
    else:
        if not np.isscalar(value):
            # sign = value / np.abs(value)
            order = np.floor(np.log10(np.nanmax(np.abs(value))) / 3 + order_bias) * 3
        else:
            order = np.floor(np.log10(np.abs(value)) / 3 + order_bias) * 3
        if order < prefix_range[0]:
            order = prefix_range[0]
        if order > prefix_range[1]:
            order = prefix_range[1]
        scaling = 10**order
        value_out = value / scaling
        prefix_out = prefixes.get(order)
        
        return(value_out, prefix_out, order)



value = -0.5366e6
value = [0.5366e6, 0.0999e-6, -0.335e-4, -0.335e6]
# value = np.array([1,30,0.4,99,801])

value_out, prefix_out, order = mprefix(value, order_bias=0.05)

if np.isscalar(value_out):
    print('{value}m = {value_out:.3f} {prefix_out}m'.format(value=value, value_out=value_out, prefix_out=prefix_out))
else:
    print(value_out, prefix_out)


#%%

# for i in range(0,1000):
#     value = 1.1**i
#     print('{:5e}:'.format(value))
#     order_bias=0.05
#     value_scaled, prefix, order = mprefix(value, order_bias=order_bias)
#     print('order_bias={:}: {:3.2f} {:} '.format(order_bias, value_scaled, prefix), end = '')

#     order_bias=0.15
#     value_scaled, prefix, order = mprefix(value, order_bias=order_bias)
#     print('    order_bias={:}: {:3.2f} {:} '.format(order_bias, value_scaled, prefix), end = '')
#     time.sleep(0.05)
    
#     order_bias=0.333
#     value_scaled, prefix, order = mprefix(value, order_bias=order_bias)
#     print('    order_bias={:}: {:3.2f} {:} '.format(order_bias, value_scaled, prefix))
#     time.sleep(0.05)
    

