__author__ = 'Sergey Tomin'

from numpy import sqrt, reshape, shape, pi, exp, zeros, array, meshgrid
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ocelot.lib.genera.src.python.convolution.convolution_gauss import convolution_1D_cpp,  convolution_2D_cpp


def plot3D_data(data, x, y):
    X,Y = meshgrid(x,y)
    fig = plt.figure()
    ax = Axes3D(fig)
    #ax = fig.add_subplot(111, projection = "3d")
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, cmap=cm.jet)

#def conditions_emitt_spread(screen):
#    if screen.ne ==1 and (screen.nx and screen.ny):
#        effect = 1
#    elif screen.ne ==1 and (screen.nx==1 and screen.ny):
#        effect = 2
#    elif screen.ne ==1 and (screen.nx and screen.ny == 1):
#        effect = 3
#    elif screen.ne >1 and (screen.nx == 1 and screen.ny == 1):
#        effect = 4
#    else:
#        effect = 0
#    return effect

def beam_sizes_on_screen(beam, screen):
    if beam.emit_x == 0 or beam.emit_y == 0:
        print ("Emittance switched off. One of emittances is zero")
        return False
    beam.sizes()
    screen.beam_size_x = sqrt(beam.sigma_x*beam.sigma_x*1e6 + (beam.sigma_xp*screen.Distance)**2) # mm
    screen.beam_size_y = sqrt(beam.sigma_y*beam.sigma_y*1e6 + (beam.sigma_yp*screen.Distance)**2) # mm

    if screen.beam_size_x == 0:
        print ("Emittance switched off. Hor. e-beam size projection on screen is zero")
        return False

    if screen.beam_size_y == 0:
        print ("Emittance switched off. Ver. e-beam size projection on screen is zero")
        return False
    print ("Hor. e-beam size projection on screen: ", screen.beam_size_x, " mm")
    print ("Ver. e-beam size projection on screen: ", screen.beam_size_y, " mm")

    #if beam.beta_x and beam.beta_y:
    #
    #    beam.sizes()
    #    #print beam.sigma_x, beam.sigma_xp
    #    #print beam.sigma_y, beam.sigma_yp
    #    #print "Distance = ", screen.Distance, " mm"
    #    screen.beam_size_x = sqrt(beam.sigma_x*beam.sigma_x*1e6 + (beam.sigma_xp*screen.Distance)**2)
    #    screen.beam_size_y = sqrt(beam.sigma_y*beam.sigma_y*1e6 + (beam.sigma_yp*screen.Distance)**2)
    #    print "Hor. e-beam size projection on screen: ", screen.beam_size_x, " mm"
    #    print "Ver. e-beam size projection on screen: ", screen.beam_size_y, " mm"
    #    return True
    #
    #elif not beam.beta_x or not beam.beta_y:
    #
    #    print "Emittance switched off: one of the beta functions is zero "
    #    return False
    #
    #else:
    #    print "Emittance switched off: beam dimensions are zero size "
    return True

def change_size(n, step, start, half_size, acc = 3):
    n_add = 0
    accuracy = acc
    if half_size != 0:
        if step == 0:
            n_add = 3*accuracy
            n = 2*n_add + 1
            start -= half_size
            step = half_size/n_add
        else:
            n_add = (int(half_size/step)+1)
            n += 2*n_add
            start -= n_add*step
    return n,step,start, n_add

def change_sizes_screen(screen, beam):
    # all sizes in mm and mrad!!!

    ## energy spread is not correct
    # screen.fund_harm_eV -> screen.start_energy
    screen.sigma_e = 2.*beam.sigma_E/beam.E*screen.start_energy
    #print "sigma_e = ", screen.sigma_e
    screen.nx_add = 0
    screen.ny_add = 0
    screen.ne_add = 0

    if beam_sizes_on_screen(beam, screen):

        if screen.ne == 1:
            screen.nx, screen.x_step, screen.x_start, screen.nx_add = change_size(screen.nx, screen.x_step, screen.x_start, screen.beam_size_x*1.) # emittance X
            screen.ny, screen.y_step, screen.y_start, screen.ny_add = change_size(screen.ny, screen.y_step, screen.y_start, screen.beam_size_y*1.) # emittance Y
            screen.ne, screen.e_step, screen.e_start, screen.ne_add = change_size(screen.ne, screen.e_step, screen.e_start, screen.sigma_e*3.) # energy spread

        elif screen.nx ==1 and screen.ny ==1:
            if screen.theta_x !=0 and screen.theta_y != 0:
                screen.nx, screen.x_step, screen.x_start, screen.nx_add = change_size(screen.nx, screen.x_step, screen.x_start, screen.theta_x*screen.Distance) # emittance X
                screen.ny, screen.y_step, screen.y_start, screen.ny_add = change_size(screen.ny, screen.y_step, screen.y_start, screen.theta_y*screen.Distance) # emittance Y
                screen.ne, screen.e_step, screen.e_start, screen.ne_add = change_size(screen.ne, screen.e_step, screen.e_start, screen.sigma_e*3.) # energy spread
            else:
                screen.nx, screen.x_step, screen.x_start, screen.nx_add = change_size(screen.nx, screen.x_step, screen.x_start, screen.beam_size_x) # emittance X
                screen.ny, screen.y_step, screen.y_start, screen.ny_add = change_size(screen.ny, screen.y_step, screen.y_start, screen.beam_size_y) # emittance Y
                screen.ne, screen.e_step, screen.e_start, screen.ne_add = change_size(screen.ne, screen.e_step, screen.e_start, screen.sigma_e*3.) # energy spread
        else:
            print ("Emittance switched off. Screen sizes are wrong. ne >1 and (nx or ny) >1")

    elif screen.nx ==1 and screen.ny ==1:
        screen.ne, screen.e_step, screen.e_start, screen.ne_add = change_size(screen.ne, screen.e_step, screen.e_start, screen.sigma_e*3.) # energy spread
    elif screen.ne == 1:
        screen.ne, screen.e_step, screen.e_start, screen.ne_add = change_size(screen.ne, screen.e_step, screen.e_start, screen.sigma_e*3.) # energy spread
    else:
        print ("1_ Emittance switched off. Screen sizes are wrong. ne >1 and (nx or ny) >1")



def emittance_on(screen):
    Pi = [] #zeros((screen.ny - 2*screen.ny_add, screen.nx - 2*screen.nx_add))
    Sigma = []# zeros((screen.ny - 2*screen.ny_add, screen.nx - 2*screen.nx_add))

    #print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    nx_ny = screen.nx*screen.ny
    #print nx_ny, shape(screen.Pi), shape(screen.Sigma), screen.beam_size_x ,screen.beam_size_y, screen.ne
    for ie in range(screen.ne):
        #print ie
        data_pi = screen.Pi[nx_ny*ie : nx_ny*(ie+1)]
        data_sigma = screen.Sigma[nx_ny*ie : nx_ny*(ie+1)]
        #print shape(data_pi), shape(data_sigma)
        data_pi = convolution_2D_cpp(data_pi,screen.Xph,screen.Yph,screen.nx_add, screen.ny_add, screen.beam_size_x ,screen.beam_size_y )
        data_sigma = convolution_2D_cpp(data_sigma,screen.Xph,screen.Yph, screen.nx_add, screen.ny_add, screen.beam_size_x ,screen.beam_size_y )

        #mean_pi = sum(data_pi)/len(data_pi)
        #data_pi[:] = mean_pi
        #mean_sigma = sum(data_sigma)/len(data_sigma)
        #data_sigma[:] = mean_sigma

        data_pi = reshape(data_pi, (screen.ny , screen.nx ))
        data_sigma = reshape(data_sigma, (screen.ny , screen.nx ))
        #print shape(data_sigma)
        data_pi = data_pi[screen.ny_add : -screen.ny_add, screen.nx_add : -screen.nx_add]
        data_sigma = data_sigma[screen.ny_add : -screen.ny_add, screen.nx_add : -screen.nx_add]
        Pi.append(data_pi)
        Sigma.append(data_sigma)

    screen.Xph = screen.Xph[screen.nx_add : screen.nx-screen.nx_add]
    screen.Yph = screen.Yph[screen.ny_add : screen.ny-screen.ny_add]
    screen.x_start += screen.nx_add*screen.x_step
    screen.y_start += screen.ny_add*screen.y_step
    screen.nx -= 2*screen.nx_add
    screen.ny -= 2*screen.ny_add
    return Pi, Sigma

def spread_on(screen, Pi,Sigma):
    sigma_e = screen.sigma_e
    ## energy spread is not correct
    # screen.fund_harm_eV -> screen.start_energy
    k = lambda energy: exp(-(energy - screen.start_energy)**2/(2.*sigma_e*sigma_e))/(sqrt(2.*pi)*sigma_e)

    if screen.ne - 2*screen.ne_add == 1:

        e_pi = zeros(shape(Pi[0]))
        e_sigma = zeros(shape(Pi[0]))

        d_e = screen.Eph[1] - screen.Eph[0]
        for ie in range(screen.ne):

            eph = screen.Eph[ie]
            #print ie, k(eph), eph
            #print array(Pi[ie])
            #print array(e_pi)
            e_pi += array(Pi[ie])*k(eph)*d_e
            e_sigma += array(Sigma[ie])*k(eph)*d_e
            #print Sigma[ie]
            #plt.plot(range(len(Sigma[ie])), Sigma[ie])
            #plt.show()
        screen.Pi = reshape(e_pi, screen.nx*screen.ny)
        #plot3D_data(e_pi, screen.Xph, screen.Yph)

        screen.Sigma = reshape(e_sigma,  screen.nx*screen.ny)
        #plot3D_data(e_sigma, screen.Xph, screen.Yph)

        screen.Total = screen.Sigma + screen.Pi
        #plot3D_data(screen.Total, screen.Xph, screen.Yph)
        #plt.show()
        #!!!!!!
        screen.Eph = screen.Eph[screen.ne_add : screen.ne - screen.ne_add]
        screen.ne -= 2*screen.ne_add
    else:
        #print shape(Pi), shape(reshape(Pi,shape(Pi)[0]))
        Pi = reshape(Pi,shape(Pi)[0])
        Pi = convolution_1D_cpp(Pi,screen.Eph, screen.ne_add, screen.sigma_e)
        screen.Pi = Pi[screen.ne_add : screen.ne - screen.ne_add]
        #print "%%%% = ",shape(Pi)
        Sigma = reshape(Sigma,shape(Sigma)[0])
        Sigma = convolution_1D_cpp(Sigma,screen.Eph, screen.ne_add, screen.sigma_e)
        screen.Sigma = Sigma[screen.ne_add : screen.ne - screen.ne_add]
        screen.Total = screen.Sigma + screen.Pi


        #print "$$$ Total= ", screen.ne, shape(screen.Total)
        screen.Eph = screen.Eph[screen.ne_add : screen.ne - screen.ne_add]
        screen.ne -= 2*screen.ne_add


def convolution_all(screen):
    #ypoint*xpoint*je + xpoint*jy + jx
    if  not screen.nx_add and not screen.ny_add and not screen.ne_add:
        print ("no convolution")
        return 0

    if screen.nx_add and screen.ny_add: #emittance switch on!!
        Pi, Sigma = emittance_on(screen)
    else:
        Pi = reshape(screen.Pi, (screen.ne, screen.nx*screen.ny))

        Sigma = reshape(screen.Sigma, (screen.ne, screen.nx*screen.ny))
    #print "shape = ", shape(Pi)
    if screen.ne_add:
        print ("energy spread on")
        spread_on(screen, Pi,Sigma)
        #plot3D_data(screen.Total, screen.Xph, screen.Yph)
        #plt.show()
        return 0
    screen.Pi = reshape(Pi,shape(Pi)[0]*shape(Pi)[1]*shape(Pi)[2])
    screen.Sigma = reshape(Sigma,shape(Sigma)[0]*shape(Sigma)[1]*shape(Sigma)[2])
    screen.Total = screen.Pi + screen.Sigma
    """
    else:
        if shape(Pi)[1] == 1 and shape(Pi)[2] == 1:
            screen.Pi = reshape(Pi,shape(Pi)[0])
            screen.Sigma = reshape(Sigma,shape(Sigma)[0])
        else:
            screen.Pi = array(Pi)
            screen.Sigma = array(Sigma)
        screen.Total = screen.Pi + screen.Sigma
    """
    #print "Total = ", shape(screen.Total)
if __name__ == "__main__":
    from radiation.em_screen import EMScreen
    from classes.screen import Screen
    screen = Screen()
    screen.fund_harm_eV = 8200
    screen.z = 300.
    screen.x = 0.00

    screen.size_x = 0.003
    screen.size_y = 0.003
    screen.nx = 1
    screen.ny = 11
    screen.start_energy = 8000
    screen.end_energy =  8400
    screen.num_energy = 1
    em_screen = EMScreen(screen)
    from classes.beam import Beam
    beam = Beam()
    beam = Beam(x=0.,xp=-0.000115*0.,y=0,yp=0)
    beam.E = 17.5
    beam.sigma_E = 0.00*17.5
    beam.Q = 0.001
    beam.emit_x = 1.752e-11
    beam.emit_y = 1.752e-10
    beam.beta_x = 33.7
    beam.beta_y = 23.218
    beam.alpha_x = 1.219
    beam.alpha_y = -0.842
    print (em_screen.nx, em_screen.x_start, em_screen.x_step)
    print (em_screen.ny, em_screen.y_start, em_screen.y_step)
    print (em_screen.ne, em_screen.e_start, em_screen.e_step)
    change_sizes_screen(em_screen, beam)
    print (em_screen.nx, em_screen.x_start, em_screen.x_step, em_screen.nx_add)
    print (em_screen.ny, em_screen.y_start, em_screen.y_step, em_screen.ny_add)
    print (em_screen.ne, em_screen.e_start, em_screen.e_step, em_screen.ne_add)
