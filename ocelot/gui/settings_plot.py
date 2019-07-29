
import matplotlib
import logging

# from pylab import rc, rcParams #tmp
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy

#in order to run decorators properly
import functools

_logger = logging.getLogger(__name__)

my_viridis = deepcopy(matplotlib.pyplot.get_cmap('viridis')) 
my_viridis.set_under('w')
def_cmap = my_viridis

def_cmap = 'viridis'
# def_cmap = 'Greys'

fntsz = 4
params = {'image.cmap': def_cmap, 'backend': 'ps', 'axes.labelsize': 3 * fntsz, 'font.size': 3 * fntsz, 'legend.fontsize': 4 * fntsz, 'xtick.labelsize': 4 * fntsz,  'ytick.labelsize': 4 * fntsz, 'text.usetex': False}
rcParams.update(params)
# plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
# rcParams["savefig.directory"] = os.chdir(os.path.dirname(__file__)) but __file__ appears to be genesis_plot
matplotlib.pyplot.ioff() #turn off interactive mode

# check if Xserver is connected
plotting_error=None
try:
    import _tkinter
    _tkinter.create()
except:
    if not "DISPLAY" in os.environ:
        plotting_error = 'Cannot plot figures: Xserver is not connected (Putty -> X11 forwarding)'
        # _logger.error('Cannot plot figures: Xserver is not connected (Putty -> X11 forwarding)')
    else:
        plotting_error = 'Cannot plot figures: Unable to connect to forwarded X server (?)'
        # _logger.error('Cannot plot figures: Unable to connect to forwarded X server (?)')

if plotting_error is not None:
    _logger.error(plotting_error)

# # re-check
# exitval = os.system('python -c "import matplotlib.pyplot as plt; plt.figure()"')
# havedisplay = (exitval == 0)
# if not havedisplay:
# # force matplotlib not ot use Xwindows backend. plots may still be plotted into e.g. *.png
# matplotlib.use('Agg')

#decorator
def if_plottable(plotting_func):
    
    def empty_func(*args, **kwargs):
        return
    
    @functools.wraps(plotting_func)
    def wrapper(*args, **kwargs):
        if plotting_error is None:
            fig = plotting_func(*args, **kwargs)
            return fig
        else:
            _logger.warning(plotting_error)
            return empty_func()

    return wrapper