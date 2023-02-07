# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: sserkez
"""

import numpy as np
import time
from ocelot.common.math_op import mprefix

value = 0.5366e6
# value = [0.5366e6, 0.0999e-6]
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
    

