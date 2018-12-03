
import matplotlib as mpl
mpl.pyplot.close('all')

set_current = "manuscript"

def current(cr='manuscript'):

    if set_current is "manuscript":
        manuscript()
    elif set_current is "presentation":
        presentation()
    else:
        default()

def manuscript():
    mpl.rcParams['toolbar'] = 'None'
    
    mpl.rcParams['font.family']     = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Helvetica' 
    mpl.rcParams['font.monospace']  = 'Courier'
    mpl.rcParams['text.color']      = 'black'
    
    mpl.rcParams['axes.facecolor']  = 'white' 
    mpl.rcParams['axes.edgecolor']  = 'black' 
    mpl.rcParams['axes.linewidth']  = 1     
    mpl.rcParams['axes.grid']       = True  
    mpl.rcParams['axes.titlesize']  = 'x-large' 
    mpl.rcParams['axes.labelsize']  = 'x-large' 
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['axes.axisbelow']  = True 
    #mpl.rcParams['axes.color_cycle'] = ['348ABD', 'A60628', '7A68A6',
    #                                    '467821', 'CF4457', '188487',
    #                                    'E24A33']

    mpl.rcParams['xtick.major.size'] = 4
    mpl.rcParams['xtick.minor.size'] = 2
    mpl.rcParams['xtick.major.pad']  = 4
    mpl.rcParams['xtick.minor.pad']  = 4
    mpl.rcParams['xtick.color']      = 'k'
    mpl.rcParams['xtick.labelsize']  = 'large'
    mpl.rcParams['xtick.direction']  = 'in'
    
    mpl.rcParams['ytick.major.size'] = 4
    mpl.rcParams['ytick.minor.size'] = 2
    mpl.rcParams['ytick.major.pad']  = 4
    mpl.rcParams['ytick.minor.pad']  = 4
    mpl.rcParams['ytick.color']      = 'k'     
    mpl.rcParams['ytick.labelsize']  = 'large'
    mpl.rcParams['ytick.direction']  = 'in'

    mpl.rcParams['grid.color']      = '0.5'
    mpl.rcParams['grid.linestyle']  = ':' 
    mpl.rcParams['grid.linewidth']  = 0.5 

    mpl.rcParams['legend.fancybox']  = False
    mpl.rcParams['legend.isaxes']    = True
    mpl.rcParams['legend.numpoints'] = 2
    mpl.rcParams['legend.fontsize']  = 'medium'
    
    mpl.rcParams['figure.figsize']  = (11, 8)   
    mpl.rcParams['figure.dpi']      = 100     
    mpl.rcParams['figure.facecolor'] = 'w' 
    mpl.rcParams['figure.edgecolor'] = 'w' 
    
    mpl.rcParams['savefig.dpi']       = '300'   
    mpl.rcParams['savefig.facecolor'] = 'white' 
    mpl.rcParams['savefig.edgecolor'] = 'white'


def presentation():
    mpl.rcParams['toolbar'] = 'None'
    
    mpl.rcParams['font.sans-serif'] = 'Arial' 
    mpl.rcParams['font.monospace']  = 'Courier'
    mpl.rcParams['text.color']      = 'white'
    mpl.rcParams['font.family']     = 'sans-serif'

    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.rm']      = 'sans:bold'
    
    mpl.rcParams['axes.facecolor']  = 'black' 
    mpl.rcParams['axes.edgecolor']  = 'white' 
    mpl.rcParams['axes.linewidth']  = 1.5     
    mpl.rcParams['axes.grid']       = True  
    mpl.rcParams['axes.titlesize']  = 'x-large' 
    mpl.rcParams['axes.labelsize']  = 'x-large' 
    mpl.rcParams['axes.labelcolor'] = 'white'
    mpl.rcParams['axes.axisbelow']  = True 
    mpl.rcParams['axes.color_cycle'] = ['348ABD', '7A68A6', 'A60628',
                                        '467821', 'CF4457', '188487',
                                        'E24A33']


    mpl.rcParams['lines.linewidth'] = 3


    mpl.rcParams['xtick.major.size'] = 4
    mpl.rcParams['xtick.minor.size'] = 2
    mpl.rcParams['xtick.major.pad']  = 4
    mpl.rcParams['xtick.minor.pad']  = 4
    mpl.rcParams['xtick.color']      = 'w'
    mpl.rcParams['xtick.labelsize']  = 'large'
    mpl.rcParams['xtick.direction']  = 'in'
    
    mpl.rcParams['ytick.major.size'] = 4
    mpl.rcParams['ytick.minor.size'] = 2
    mpl.rcParams['ytick.major.pad']  = 4
    mpl.rcParams['ytick.minor.pad']  = 4
    mpl.rcParams['ytick.color']      = 'w'     
    mpl.rcParams['ytick.labelsize']  = 'large'
    mpl.rcParams['ytick.direction']  = 'in'

    mpl.rcParams['grid.color']      = '0.5'
    mpl.rcParams['grid.linestyle']  = '-' 
    mpl.rcParams['grid.linewidth']  = 1 

    mpl.rcParams['legend.fancybox']  = False
    mpl.rcParams['legend.isaxes']    = True
    mpl.rcParams['legend.numpoints'] = 2
    mpl.rcParams['legend.fontsize']  = 'medium'
    
    mpl.rcParams['figure.figsize']  = (11, 8)   
    mpl.rcParams['figure.dpi']      = 100     
    mpl.rcParams['figure.facecolor'] = 'k' 
    mpl.rcParams['figure.edgecolor'] = 'k' 
    
    mpl.rcParams['savefig.dpi']       = '200'   
    mpl.rcParams['savefig.facecolor'] = 'black' 
    mpl.rcParams['savefig.edgecolor'] = 'black'

