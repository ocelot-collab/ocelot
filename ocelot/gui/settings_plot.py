
import matplotlib, logging, os

# from pylab import rc, rcParams #tmp
from matplotlib import rc, rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy
from ocelot.common.ocelog import ind_str

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

#matplotlib.pyplot.style.use("dark_background") White axes on dark background
# matplotlib.use('Agg')

# check if Xserver is connected
plotting_error=None
#try:
#    import _tkinter
#    _tkinter.create()
#except:
#    if not "DISPLAY" in os.environ:
#        plotting_error = 'Cannot plot figures: Xserver is not connected (Putty -> X11 forwarding)'
#        # _logger.error('Cannot plot figures: Xserver is not connected (Putty -> X11 forwarding)')
#    else:
#        plotting_error = 'Cannot plot figures: Unable to connect to forwarded X server (?)'
#        # _logger.error('Cannot plot figures: Unable to connect to forwarded X server (?)')

#if plotting_error is not None:
#    _logger.error(plotting_error)

# # re-check
# exitval = os.system('python -c "import matplotlib.pyplot as plt; plt.figure()"')
# havedisplay = (exitval == 0)
# if not havedisplay:
# # force matplotlib not ot use Xwindows backend. plots may still be plotted into e.g. *.png
# matplotlib.use('Agg')


#decorator
def save_figure(plotting_func):
    
    @functools.wraps(plotting_func)
    def wrapper(*args, **kwargs):
        savepath = kwargs.pop(savepath, None)
        fig = plotting_func(*args, **kwargs)
        if savepath is not None:
            fig.savefig(savepath, format=savepath.split('.')[-1])
        return fig
    
    return wrapper


def save_show(plotting_func):
    @functools.wraps(plotting_func)
    def wrapper(*args, **kwargs):
        savefig = kwargs.pop('savefig', False)
        showfig = kwargs.pop('showfig', True)
        closefig = kwargs.pop('closefig', True)
        fig = plotting_func(*args, **kwargs)
        
        matplotlib.pyplot.draw()
        if savefig != False:
            if savefig == True:
                savefig = 'png'
            if savefig in ['png' , 'eps', 'pdf', 'jpeg']:
                savepath = edist.filePath + '.' + savefig
            else:
                savepath = savefig
            _logger.debug(ind_str + 'saving to {}'.format(savepath))
            matplotlib.pyplot.savefig(savepath, format=savepath.split('.')[-1])
        
        if showfig:
            # rcParams["savefig.directory"] = os.path.dirname(edist.filePath)
            matplotlib.pyplot.show()
        if closefig:
            # # plt.close('all')
            matplotlib.pyplot.close(fig)
        
        # if savepath is not None:
            # fig.savefig(savepath, format=savepath.split('.')[-1])
        return fig
    return wrapper


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