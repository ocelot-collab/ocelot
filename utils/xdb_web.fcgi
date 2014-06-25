#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# /usr/lib/cgi-bin/test.fcgi
# sudo service apache2 restart

from cgi import escape
import sys, os
from flup.server.fcgi import WSGIServer

import xframework.adaptors.genesis as genesis
from xframework.utils.xdb import Xdb
from xframework.cpbd.elements import Magnet, Quadrupole, RBend, Drift, Undulator, MagneticLattice
from xframework.cpbd.beam import Beam
from xframework.cpbd.optics import *

import numpy as np
from numpy import *

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import cStringIO


def header():
    
    ss = '<script type="text/javascript"> var userAgent = navigator.userAgent; document.write(userAgent); </script>'
    css = '<link rel="stylesheet" type="text/css" href="/css/xdb.css" />'
    
    return '<!DOCTYPE html><html><head><title>XDB</title>'+css+'</head><body><div><a href="/cgi-bin/test.fcgi/browse">Browse</a><br><br></div>' 

def footer():
    return '<!--<tiny> (C)Ilya.Agapov@xfel.eu 2013 </tiny>--></body></html>' 


def undulators(xdb):
    
    try:
        s = '<h4>Undulators</h4><table border="0">'
        g = xdb.file['Undulators']
        for u in g.keys():
            s = s +  '<tr><td>' + str(u) + '</td><td><a href = "show/input/'+str(u)+'">input</a></td><td><a href = "show/twiss/'+str(u)+'">twiss</a></td>'

            for u2 in g[u].keys():
                s = s + '<td>' +  str(u2) + '</td>'
                
            s += '</tr>'
        return str(s) + '</table>'
    except:
        return 'error retrieving undulator list'

def beams(xdb):
    
    try:
        global s
        s = '<h4>Beams</h4><table border="0">'
        global g
        g = xdb.file['Beams']
        def printname(name):
            global s, g
            if len( g[name].keys()) >= 1:
                s = s + '<tr><td>' +  name + '</td></tr>\n'
            else:
                s = s + '<tr><td><a href="http://www.desy.de">' +  name + '</a></td></tr>'
        g.visit(printname)
                
        return str(s) + '</table>'
    except:
        return 'error retrieving undulator list'


def fel_parameters(xdb):
    try:
        global s
        s = '<h4>FEL parameters::</h4><table border="0">'
        g = xdb.file['FEL']
        
        for u in g.keys():
            g2 = g[u]
            for u2 in g2.keys():
                g3 = g2[u2]
                for u3 in g3.keys():
                    g4 = g3[u3]
                    for u4 in g4.keys():
                        g5 = g4[u4]
                        for u5 in g5.keys():
                            name = str(u) + '/' + str(u2) + '/' + str(u3) + '/' + str(u4) + '/' + str(u5)
                            s = s +  '<tr><td><a href="/cgi-bin/test.fcgi/show_stat/'+name+'/v1">' + name + '</a></td></tr>'
             
        return str(s) + '</table>'
    except:
        return 'error retrieving fel parameter list'
    
def input(xdb, undulator):
    try:
        s = xdb.read_undulator_config(undulator)
        return '<pre>' + str(s) + '</pre>'
    except:
        return 'error retrieving input for: ' + undulator 

def app2(environ, start_response):
    start_response('200 OK', [('Content-Type', 'text/html')])

    yield '<h1>FastCGI Environment</h1>\n\n'
    yield '<table>'
    for k, v in sorted(environ.items()):
         yield '<tr><th>%s</th><td>%s</td></tr>' % (escape(k), escape(v))
    yield '</table>'
    #yield str(sys.path)

def parse_path(path):
    path = path.split('/')[1:]
    parsed_path = {}
    parsed_path['action'] = path[0]
    
    if parsed_path['action'] in ('img','show'):
        parsed_path['data_type'] = path[1]
        
        if parsed_path['data_type'] == 'input':
            undulator = path[2]
            parsed_path['undulator'] = undulator
        else:
            undulator = ''
            for u in path[2:]:
                undulator = undulator + u + '/'
            parsed_path['undulator'] = undulator

    return parsed_path
    
    
def app(environ, start_response):
    
    #start_response('200 OK', [('Content-Type', 'text/html')])
    
    xdb = Xdb(index_file='/home/iagapov/data/xdb/index.h5', mode='r')
    
    path = environ['PATH_INFO']
    
    p_path = parse_path(path)
                
    if p_path['action'] == 'browse': 
        start_response('200 OK', [('Content-Type', 'text/html')])
        yield header()
        yield undulators(xdb)
        yield fel_parameters(xdb)
        #yield beams(xdb)
        #yield correlation_plots(xdb)
        yield footer()
        
        
    if p_path['action'] == 'show':
        
        if p_path['data_type'] == 'input':
            start_response('200 OK', [('Content-Type', 'text/html')])
            yield header()
            yield input(xdb,p_path['undulator'])
            yield footer()
        else:
            start_response('200 OK', [('Content-Type', 'text/html')])
            yield header()
            path = path.replace('/show/','/img/')
            yield '<img src="/cgi-bin/test.fcgi'+path+'"></img>'
            yield footer()

    if p_path['action'] == 'show_stat':
        
        start_response('200 OK', [('Content-Type', 'text/html')])
        yield header()
        path1 = path.replace('/show_stat/','/img/power_z/') 
        yield '<img src="/cgi-bin/test.fcgi'+path1+'"></img><br>'
        path2 = path.replace('/show_stat/','/img/pulse/') 
        yield '<img src="/cgi-bin/test.fcgi'+path2+'"></img><br>'
        path3 = path.replace('/show_stat/','/img/spec/') 
        yield '<img src="/cgi-bin/test.fcgi'+path3+'"></img><br>'        
        yield footer()

        
    if p_path['action'] == 'img': 
        
        start_response('200 OK', [('Content-Type', 'image/png')])
        fig=Figure()
        ax=fig.add_subplot(111)
        ax.grid(True)

        if p_path['data_type'] == 'twiss':
                f = xdb.read_undulator_config(p_path['undulator'])
    
                exec(f)
                
                lat = MagneticLattice(sase)
    
                tw0 = Twiss(beam)
                print tw0
                tws=twiss(lat, tw0, nPoints = 1000)
                p1, = ax.plot(map(lambda p: p.s, tws), map(lambda p: p.beta_x, tws), lw=2.0)
                p2, = ax.plot(map(lambda p: p.s, tws), map(lambda p: p.beta_y, tws), lw=2.0)
                ax.legend([p1,p2], [r'$\beta_x$',r'$\beta_y$'])

        
        if p_path['data_type'] == 'power_z':
                                    
            p_med = xdb.file['FEL/' + p_path['undulator'] + '/power_z_med']
            p_mean = xdb.file['FEL/' + p_path['undulator'] + '/power_z_mean']
            p_std = xdb.file['FEL/' + p_path['undulator'] + '/power_z_std']
            try:
                z = xdb.file['FEL/' + p_path['undulator'] + '/z']
            except:
                z = np.arange(0,len(p_med))
                
            ax.plot(z, p_med, '-', lw=3)
            ax.errorbar(z, p_mean, yerr=p_std, fmt='g--',lw=1, capsize=5)

           
        if p_path['data_type'] == 'pulse':
           
            p_med = xdb.file['FEL/' + p_path['undulator'] + '/pulse_med']
            p_mean = xdb.file['FEL/' + p_path['undulator'] + '/pulse_mean']
            p_std = xdb.file['FEL/' + p_path['undulator'] + '/pulse_std']
            t = xdb.file['FEL/' + p_path['undulator'] + '/t']
            ax.plot(t, p_med, '-', lw=3)
            ax.errorbar(t, p_mean, yerr=p_std, fmt='g--',lw=1, capsize=5)

        if p_path['data_type'] == 'spec':
           
            s_med = xdb.file['FEL/' + p_path['undulator'] + '/spec_med']
            s_mean = xdb.file['FEL/' + p_path['undulator'] + '/spec_mean']
            s_std = xdb.file['FEL/' + p_path['undulator'] + '/spec_std']
            try:
                freq_ev = xdb.file['FEL/' + p_path['undulator'] + '/freq_ev']
            except:
                freq_ev = 1000 * np.arange(0,len(s_mean))
                
            ax.plot(freq_ev, s_med, '-', lw=3)
            ax.errorbar(freq_ev, s_mean, yerr=s_std, fmt='g--',lw=1, capsize=5)



        canvas=FigureCanvas(fig)
        buf = cStringIO.StringIO()
        canvas.print_png(buf)
        yield buf.getvalue()    



WSGIServer(app).run()
