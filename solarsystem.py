import numpy as np
import sys
from numpy.linalg import norm
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem

seed = 75041 #daniehei001

system = AST1100SolarSystem(seed,hasMoons=False)


#Variables for simulations
simulationTime = 60#In years
steps = 20000 #Steps per year
dt = 1./steps

save = False

G = 4*np.pi**2
c = 63239.7263 #AU/year
AU_per_year_to_SI = 4740.57172
au_to_m = 1.496e11

use_many_body = False
number_for_N_body = 3



writeingFreq = 100

if (float(steps*simulationTime)/writeingFreq - steps*simulationTime/writeingFreq) != 0.0:
    print "change writeingFreq"
    sys.exit()


time = np.linspace(0,simulationTime,int(steps*simulationTime/writeingFreq))


print "The habitable zone is within [%g,%g] [km]" %((system.starRadius/2)*(system.temperature/260.)**2,(system.starRadius/2)*(system.temperature/390.)**2)




if (not use_many_body and (number_for_N_body != system.numberOfPlanets)):

    N = system.numberOfPlanets

    try:
        postions = np.zeros([2,N + 1,steps*simulationTime/writeingFreq]) #Number of datapoints, numb of planets and star, 2 coordinates
    except:
        print "change writeingFreq"

    temp_pos = np.zeros([2,N+1])
    velocity = np.zeros([2,N + 1])#,steps*simulationTime])
    sun_vel = np.zeros([2,steps*simulationTime/writeingFreq])
    sun_pos = np.zeros([2,steps*simulationTime/writeingFreq])

    masses = np.zeros(N+1)
    periods = np.zeros(N+1)
    acceleration = np.zeros([2,N+1])



    sunMass = system.starMass
    masses[N] = sunMass

    b = 2900 # Wien's displacement constant in micrometers
    sun_waveLenght = b/system.temperature

    print "mass: ", sunMass
    print "Temp: ",system.temperature
    print "Wavelenght: ", sun_waveLenght
    print "Radius: ",system.starRadius
    print "___"
    print len(system.mass)




    #Initializing the arrays
    for i in range(N):
        masses[i] = system.mass[i]
        print "Planet %s " %str(i+1)
        print "Planet mass: ",masses[i]
        print "Planet radius: ",system.radius[i]
        print "Planet e: ", system.e[i]
        print "Planet Atmosphere: ", system.rho0[i]
        print "Planet a[km]: ", system.a[i]*au_to_m/1000


        postions[0,i,0] = system.x0[i]
        postions[1,i,0] = system.y0[i]
        temp_pos[0,i] = system.x0[i]
        temp_pos[1,i] = system.y0[i]



        velocity[0,i] = system.vx0[i]
        velocity[1,i] = system.vy0[i]

        periods[i] = (2*np.pi/system.period[i])*360.
        print "Planet period: ",system.period[i]

        print "----"
else:

    N = number_for_N_body

    b = 2900 # Wien's displacement constant in micrometers
    sun_waveLenght = b/system.temperature
    sunMass = system.starMass

    print "mass: ", sunMass
    print "Temp: ",system.temperature
    print "Wavelenght: ", sun_waveLenght
    print "Radius: ",system.starRadius
    print "___"
    print len(system.mass)


    try:
        postions = np.zeros([2,N + 1,steps*simulationTime/writeingFreq]) #Number of datapoints, numb of planets and star, 2 coordinates
    except:
        print "change writeingFreq"

    index_biggest_masses = []
    temp_pos = np.zeros([2,N+1])
    velocity = np.zeros([2,N + 1])#,steps*simulationTime])
    sun_vel = np.zeros([2,steps*simulationTime/writeingFreq])
    sun_pos = np.zeros([2,steps*simulationTime/writeingFreq])

    complet_list_masses = np.array(system.mass)
    masses = np.zeros(N+1)
    periods = np.zeros(N+1)
    acceleration = np.zeros([2,N+1])


    masses[N] = sunMass

    index_biggest = np.argsort(complet_list_masses)
    for i in range(N):
        index_biggest_masses.append(index_biggest[-i-1])
        masses[i] = system.mass[index_biggest[-i-1]]


    for i in range(N):
        masses[i] = system.mass[index_biggest_masses[i]]
        print "Planet %s " %str(i+1)
        print "Planet mass: ",masses[i]
        print "Planet radius: ",system.radius[index_biggest_masses[i]]
        print "Planet e: ", system.e[i]
        print "Planet a: ", system.a[i]

        postions[0,i,0] = system.x0[index_biggest_masses[i]]
        postions[1,i,0] = system.y0[index_biggest_masses[i]]
        temp_pos[0,i] = system.x0[index_biggest_masses[i]]
        temp_pos[1,i] = system.y0[index_biggest_masses[i]]



        velocity[0,i] = system.vx0[index_biggest_masses[i]]
        velocity[1,i] = system.vy0[index_biggest_masses[i]]

        periods[i] = (2*np.pi/system.period[index_biggest_masses[i]])*360.
        print "Planet dtheta: ",periods[i]

        print "----"





if(use_many_body):

    mom = np.zeros(2)
    mom[0] = np.sum(velocity[0,:]*masses)
    mom[1] = np.sum(velocity[1,:]*masses)

    velocity[:,-1] = -mom/sunMass



    for i in xrange(N+1):
        acceleration[:,i] = 0
        rx = (temp_pos[0,:] - float(temp_pos[0,i]))**2
        ry = (temp_pos[1,:] - float(temp_pos[1,i]))**2
        other_p = np.power(rx+ry,3./2) > 1.e-5
        some = np.maximum(np.power(rx+ry,3./2),1.e-5)
        acceleration[0,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[0,:] - float(temp_pos[0,i]))))
        acceleration[1,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[1,:] - float(temp_pos[1,i]))))



    velocity += 0.5*acceleration*dt

    for step in xrange(steps*simulationTime-1):

        temp_pos += velocity*dt

        for i in xrange(N+1):
            acceleration[:,i] = 0
            rx = (temp_pos[0,:] - float(temp_pos[0,i]))**2
            ry = (temp_pos[1,:] - float(temp_pos[1,i]))**2
            other_p = np.power(rx+ry,3./2) > 1.e-5
            some = np.maximum(np.power(rx+ry,3./2),1.e-5)
            acceleration[0,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[0,:] - float(temp_pos[0,i]))))
            acceleration[1,i] = np.sum(other_p*(G*masses[:]/some*(temp_pos[1,:] - float(temp_pos[1,i]))))


        velocity += acceleration*dt

        if ((step+1)%writeingFreq == 0):
            postions[:,:,step/writeingFreq + 1] = temp_pos[:,:]
            sun_pos[:,step/writeingFreq +1] = temp_pos[:,-1]
            sun_vel[:,step/writeingFreq +1] = velocity[:,-1]

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""



else:
    acceleration = -G*sunMass/(norm(temp_pos[:,:N],axis = 0)**3)*temp_pos[:,:N]
    velocity[:,:N] += 0.5*acceleration[:,:N]*dt


    for step in xrange(steps*simulationTime-1):


        temp_pos[:,:N] += velocity[:,:N]*dt

        acceleration = -G*sunMass/(norm(temp_pos[:,:N],axis = 0)**3)*temp_pos[:,:N]


        velocity[:,:N] += acceleration[:,:N]*dt


        if ((step+1)%writeingFreq == 0):
            postions[:,:,step/writeingFreq + 1] = temp_pos[:,:]

        print (float(step)/(steps*simulationTime))*100, "%            \r",
    print ""



#system.checkPlanetPositions(postions,simulationTime,steps/writeingFreq)
#system.orbitXml(postions[:,:N,:],time)


plt.xlabel("Position [AU]")
plt.ylabel("Position [AU]")
plt.axis("equal")
for p in range(N):
    plt.plot(postions[0,p,-1],postions[1,p,-1],"o")
    plt.plot(postions[0,p,0],postions[1,p,0],"x")
    plt.plot(postions[0,p,:],postions[1,p,:])

plt.plot(postions[0,N,-1],postions[1,N,-1],"y*")

plt.plot(postions[0,N,:],postions[1,N,:])



plt.show()
if(save and not use_many_body):
    outFilePos = open("filePositions.npy", "wb")
    np.save(outFilePos, postions[:,:N,:])
    outFilePos.close()
if (use_many_body):
    random_noice_vel = []
    for i in sun_vel[0]:
        random_noice_vel.append(i + np.random.normal(0,(1./5)*np.max(sun_vel[0])))
    plt.xlabel("Time[years]")
    plt.ylabel("Velocity [m/s]")
    plt.plot(time,sun_vel[0]*AU_per_year_to_SI)
    plt.show()
    plt.xlabel("Time[years]")
    plt.ylabel("Velocity [m/s]")
    plt.plot(time,np.array(random_noice_vel)*AU_per_year_to_SI)
    plt.show()
