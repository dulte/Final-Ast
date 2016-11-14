# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:30:05 2016
@author: Daniel
"""
from AST1100SolarSystem import AST1100SolarSystem
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import random
import numpy as np
import matplotlib.pyplot as plt
import os

class engine:
    def __init__(self, sideSize, holeSize, fuel, temp,interval,steps,number_of_particles):
        self.L = float(sideSize)
        self.hole = holeSize
        self.fuel = fuel
        self.T = temp
        self.k = 1.38e-23
        self.dt = interval/steps
        self.steps = steps
        self.interval = interval

        self.launcher_mass = 1100

        seed = 75041
        self.myStarSystem = AST1100SolarSystem(seed)



        self.m =  3.3e-27
        self.N = number_of_particles

        self.sigma = np.sqrt(self.k*self.T/self.m)

        self.force_on_top = 0
        self.momentum_gaines = 0
        self.force_gained = 0

        self.numb_escaping = 0
        self.numb_coll = 0
        self.mass_escaped_per_sec = 0

        self.engine_stated = False
        self.G = 6.674e-11

        #self.v_escape = (2.04+1.2)*1.1*4743.7173611#4.64893633785*4743.7173611#12473

        self.solar_to_kg = 1.988435e30
        self.v_escape = np.sqrt(2*self.G*self.myStarSystem.mass[0]*self.solar_to_kg/(self.myStarSystem.radius[0]*1000.))


        self.total_time = 0

        self.fuel_list = []
        self.velocity_list = []





        self.x = np.zeros([int(self.N),3])
        self.v = np.zeros([int(self.N),3])




        for i in range(int(self.N)):
            for j in range(3):
                self.x[i,j] = -self.L/2+(random.random()*self.L)
                self.v[i,j] = random.gauss(0,self.sigma)



    def integrate(self):
        self.x += self.v*self.dt


    def detect_collision_vec(self):

        for i in range(3):

            if i == 1:
                index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2, np.logical_and(np.logical_and(self.x[:,i+1] < self.hole/2, self.x[:,i+1] > -self.hole/2),np.logical_and(self.x[:,i-1] < self.hole/2,self.x[:,i-1] > -self.hole/2)))
                #index_part_through_floor = np.logical_and(self.x[:,i] < -self.L/2., np.logical_and(self.x[:,i-1] < self.hole/2.,self.x[:,i-1] > -self.hole/2.),np.logical_and(self.x[:,i+1] < self.hole/2.,self.x[:,i+1] > -self.hole/2.))
                self.momentum_gaines += np.abs(np.sum(self.v[index_part_through_floor,i])*self.m)
                self.numb_escaping += np.sum(index_part_through_floor)
                self.force_on_top += 2*np.sum(self.v[self.x[:,i] > self.L/2.,i])*self.m/self.dt
                self.numb_coll += np.sum(self.x[:,i] < -self.L/2.)





            self.v[np.abs(self.x[:,i]) > self.L/2.,i] = -self.v[np.abs(self.x[:,i]) > self.L/2.,i]


            self.x[self.x[:,i] > self.L/2.,i] = self.L/2.
            self.x[self.x[:,i] < -self.L/2.,i] = -self.L/2.




    def print_ke(self):
        print "Calculated energy: ", (3./2)*self.k*self.T
        ke = 0
        mid = 0
        for i in range(int(self.N)):
            mid = 0
            for j in range(3):
                mid += self.v[i,j]**2

            ke += 0.5*self.m*mid

        print "Numerical energy: ", ke/self.N

    def print_press(self):
        print "Numerical Pressure: " , (self.force_on_top/(self.L**2))/steps
        print "Analytical Pressure: ", self.N*self.k*self.T/(self.L**3)

    def startEngine(self):

        for i in range(self.steps):

            self.detect_collision_vec()

            self.integrate()
            print (float(i)/self.steps)*100, "%            \r",
        print ""

        self.force_gained = self.momentum_gaines/(self.dt*self.steps)
        self.mass_escaped_per_sec = self.numb_escaping*self.m/self.interval

        self.engine_stated = True
        self.numb_escaping_per_sec = self.numb_escaping/self.interval

        print "Engine has stated!"
        print "Data from the engine: "
        print "--------------------------------"
        print "Force: ", self.force_gained
        self.print_press()
        self.print_ke()
        print "Number of particles colliding with floor: ", self.numb_coll
        print "Number of particles escaped: ",self.numb_escaping
        print "Number of particles escaped per sec: ",self.numb_escaping_per_sec
        print "Mass lost per sec: " ,self.mass_escaped_per_sec
        print "--------------------------------"


    def plot_box(self,frames):
        def update_lines(num, dataLines, lines) :
            for line, data in zip(lines, dataLines) :
                line.set_data(data[0:2, num-1:num])
                line.set_3d_properties(data[2,num-1:num])
            return lines


        # Attach 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        m = frames   # number of frames in the animation
        n = self.N   # number of particles you want to animate
        N = self.steps # number of time steps in your data

        data = np.zeros([n, 3, N]) # this is the positions of the particles

        for i in range(N):

            for p in range(self.N):
                data[p,:,i] = self.x[p,:]
            self.detect_collision_vec()
            self.integrate()
            print (float(i)/self.steps)*100, "%            \r",
        print ""


        lines = [i for i in range(n)]
        for i in range(n):
            lines[i] = [ax.plot(data[i][0,0:1],
            data[i][1,0:1], data[i][2,0:1], 'o')[0]]


        ax.set_xlim3d([-self.L/2., self.L/2.])
        ax.set_xlabel('X')

        ax.set_ylim3d([-self.L/2., self.L/2.])
        ax.set_ylabel('Y')

        ax.set_zlim3d([-self.L/2., self.L/2.])
        ax.set_zlabel('Z')

        ax.set_title('3D random particles')


        ani = [i for i in range(n)]
        for i in range(n):
            ani[i] = animation.FuncAnimation(fig,
            update_lines, m, fargs=([data[i]], lines[i]),
            interval=50, blit=False)
        plt.show()




    def launch_simpel(self, numberOfBoxes, burnTime, launch_step):

        launch_dt = float(burnTime)/launch_step

        self.total_mass_escaping_per_sec = numberOfBoxes*self.mass_escaped_per_sec
        mass_fuel = numberOfBoxes*self.N*self.m
        time = 0
        self.total_force = numberOfBoxes*self.force_gained
        fuel_mass = []
        fuel_mass.append(0)
        acceleration = 0

        reached_escape_vel = False


        launcher_vel = []
        launcher_vel.append(0)

        while (time < burnTime and launcher_vel[-1] < self.v_escape):

            acceleration = self.total_force/(self.launcher_mass)# + mass_fuel)

            launcher_vel.append(launcher_vel[-1] + acceleration*launch_dt)

            fuel_mass.append(fuel_mass[-1] + self.total_mass_escaping_per_sec*launch_dt)

            time += launch_dt

            print (float(time)/burnTime)*100, "%            \r",
        print ""

        if (launcher_vel[-1] > self.v_escape):
            reached_escape_vel = True

        self.checkFuelAmount(numberOfBoxes,fuel_mass[-1])

        if(reached_escape_vel):
            print "Fuel needed: ", fuel_mass[-1]
            print self.mass_escaped_per_sec/self.force_gained
            print "Analytical: ", self.total_mass_escaping_per_sec*1200/(1-np.exp(self.force_gained/(self.v_escape*self.mass_escaped_per_sec))) - self.launcher_mass
            print "After ", time, " seconds!"
        else:
            print "Did not reach escape velocity"
            print "Only reach ", launcher_vel[-1], " m/s"


        x = np.linspace(0,time, len(launcher_vel))*launch_dt*60
        plt.xlabel("Time[s]")
        plt.ylabel("Fuel[kg],Velocity[m/s]")
        plt.plot(x,[self.v_escape for i in range(len(launcher_vel))])
        plt.plot(x,launcher_vel)
        plt.plot(x,fuel_mass)
        plt.legend(["Escape Velocity","Velocity","Mass of fuel"], loc = 2)
        plt.show()


    def launch(self, numberOfBoxes, burnTime, launch_step):

        launch_dt = float(burnTime)/launch_step

        self.total_mass_escaping_per_sec = numberOfBoxes*self.mass_escaped_per_sec
        time = 0
        self.total_force = numberOfBoxes*self.force_gained
        fuel_mass = []
        fuel_mass.append(self.fuel)
        acceleration = 0

        self.reached_escape_vel = False


        launcher_vel = []
        launcher_vel.append(0)

        while (time < burnTime and launcher_vel[-1] < self.v_escape):

            acceleration = self.total_force/(self.launcher_mass + fuel_mass[-1])

            launcher_vel.append(launcher_vel[-1] + acceleration*launch_dt)

            fuel_mass.append(fuel_mass[-1] - self.total_mass_escaping_per_sec*launch_dt)

            time += launch_dt

            print (float(time)/burnTime)*100, "%            \r",
        print ""

        if (launcher_vel[-1] > self.v_escape):
            self.reached_escape_vel = True

        if(self.reached_escape_vel):
            print "Reached escape velocity with: ", fuel_mass[-1], " to spare!"
            #print "Fuel needed: ", fuel_mass[-1]
            print "After ", time, " seconds!"
            print "Analytical velocity: ",(self.force_gained/self.mass_escaped_per_sec)*np.log((self.launcher_mass+self.fuel)/((self.launcher_mass+self.fuel) - self.total_mass_escaping_per_sec*time))
            print "Analytical fuel used: ", np.exp(self.mass_escaped_per_sec*self.v_escape/self.force_gained)*self.total_mass_escaping_per_sec*time/(np.exp(self.v_escape*self.mass_escaped_per_sec/self.force_gained)-1) - self.launcher_mass
            print "Used: ", self.fuel - fuel_mass[-1]
            self.fuel = fuel_mass[-1]
            self.fuel_list = fuel_mass
            self.velocity_list = launcher_vel
            self.total_time += time

        else:
            print "Did not reach escape velocity, fell down and crashed!"
            print "Only reach ", launcher_vel[-1], " m/s"
            return

        x = np.linspace(0,time, len(launcher_vel))*launch_dt*60
        plt. xlabel("Time[s]")
        plt.ylabel("Fuel[kg],Velocity[m/s]")
        plt.plot(x,[self.v_escape for i in range(len(launcher_vel))])
        plt.plot(x,launcher_vel)
        plt.plot(x,fuel_mass)
        plt.legend(["Escape Velocity","Velocity","Mass of fuel"], loc = 2)
        plt.show()


    def checkFuelAmount(self, numberOfBoxes,fuel):


        self.myStarSystem.massNeededCheck(numberOfBoxes,self.v_escape,self.force_gained, self.numb_escaping/interval, fuel)



    def boost_fuel_check(self, boost, current_mass_of_fuel):
        if (not self.reached_escape_vel):
            print "It seems you havnt lauched yet"
            return


        dt = 1./10000
        time = 0
        vel = 0

        while (vel < boost and current_mass_of_fuel > 0):
            vel += (self.total_force/(self.launcher_mass+current_mass_of_fuel))*dt
            current_mass_of_fuel -= self.total_mass_escaping_per_sec*dt
            time += dt
            if (current_mass_of_fuel < 0):
                print "he"

        if vel < boost:
            print "It seem you didnt have enough fuel"
            print "You only got ", vel, "m/s"
        else:
            print "You reached your desired velocity with ", current_mass_of_fuel, " fuel left"
            print "After ", time, " sec"



    def boost(self, boost):
        if (not self.reached_escape_vel):
            print "It seems you havnt lauched yet"
            return
        if (self.fuel_list[-1] == 0):
            print "You are out of fuel, and are stranded in space..."

        dt = 1./100000
        time = 0


        while (vel < boost and current_mass_of_fuel > 0):
            self.velocity_list( self.velocity_list[-1] + self.total_force/(self.launcher_mass+self.fuel_list[-1]))
            self.fuel_list.append(self.fuel_list[-1] - self.total_mass_escaping_per_sec*dt)
            time += dt

        if vel < boost:
            print "It seem you didnt have enough fuel"
            print "You only got ", self.velocity_list[-1], "m/s"
        else:
            print "You reached your desired velocity with ", self.fuel_list[-1], " fuel left"
            print "After ", t, " sec"



        self.total_time += time


    def check_data(self):
        print "Current velocity: ", self.velocity_list[-1], " m/s"
        print "Current fuel: ", self.fuel_list[-1]
        print "Time into mission: ", self.total_time


    def calcNumbOfBoxes(self,time,force = 0):
        if force == 0:
            return self.v_escape*self.launcher_mass/(self.force_gained*time)
        else:
            return self.v_escape*self.launcher_mass/(force*time)





interval = 1e-9
steps = 1000
frames = 1000
L = 1e-6
H = L/2
T = 10000
numb_part = 1e5 #25

e = engine(L,H,4100,T,interval,steps,numb_part)


#e.plot_box(frames)

e.startEngine()
e.launch(1.333e13, 40*60,100000)
e.launch_simpel(1.333e13, 40*60,100000)


print "Number of boxes: ", e.calcNumbOfBoxes(10*60)


fuel = e.launcher_mass*(np.exp(e.mass_escaped_per_sec*e.v_escape/e.force_gained) - 1)#np.exp(e.mass_escaped_per_sec*e.v_escape/e.force_gained)*e.total_mass_escaping_per_sec*413/(np.exp(e.v_escape*e.mass_escaped_per_sec/e.force_gained)-1) - e.launcher_mass
print "Fuel needed", fuel
e.boost_fuel_check(e.v_escape,fuel)
