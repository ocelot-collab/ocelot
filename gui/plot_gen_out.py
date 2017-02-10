
import os.path
from ocelot.gui.genesis_plot import *
from ocelot.utils.xfel_utils import *

# dirrectory with results (with run_XXX folders)
results_dir = '/data/netapp/xfel/yevgeniy/moga_test/moga_fit_new_results/iter_0/'

i = 0
dir = results_dir + 'run_' + str(i) + '/'

while(os.path.exists(dir)):
    
    file = dir + 'run.' + str(i) + '.s1.gout'
    if(file): 
        background('''plot_gen_out_all("'''+file+'''", savefig='png', choice=(1,1,1,1,6.14,0,0,0,0,0,0), showfig=True)''')
    
    i += 1
    dir = results_dir + 'run_' + str(i) + '/'

plot_gen_stat(proj_dir=results_dir, run_inp=[], stage_inp=[], param_inp=[], s_param_inp=['p_int','energy','r_size_weighted'], z_param_inp=[], dfl_param_inp=[], s_inp=['max'], z_inp=[0,'end'], savefig=1, saveval=1, showfig=0, debug=0)
