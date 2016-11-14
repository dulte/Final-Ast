import numpy as np
import matplotlib.pyplot as plt
from AST1100SolarSystem import AST1100SolarSystem
import scipy.interpolate as inter
from numpy.linalg import norm
import time as tid
import seaborn
import sys


class satelite:
    def __init__(self,destination):
        interval = 1e-9
        steps = 1000
        frames = 1000



        self.G = 4*np.pi**2

        self.km_to_au = 6.684587e-9



        #self.seed = 69558
        self.seed = 75041


        self.system = AST1100SolarSystem(self.seed)



        self.velocity = np.zeros(2)
        self.position = np.zeros(2)
        self.sim_time = 2.15
        self.steps_per_year = 30000#365*24*60#300000
        self.dt = 1./self.steps_per_year




        self.writeingFreq = 10

        self.pos_over_time = np.zeros ((2,int(self.sim_time*self.steps_per_year/self.writeingFreq)))#np.zeros((2,self.sim_time*self.steps_per_year/self.writeingFreq))

        self.time_planet_simulated = 20

        self.use_time0 = False


        self.save = False








        self.numberOfPlanets = self.system.numberOfPlanets
        self.destination_planet = destination
        self.home_planet = 0

        self.escape_velocity = self.calc_escape()



        self.planetMasses = np.array(self.system.mass)
        self.starMass = self.system.starMass
        self.starRadius = self.system.starRadius
        self.starTemp = self.system.temperature


        with open("planetPositions.npy", "rb") as npy:
            [self.planetPositions,self.time] = np.load(npy)



        self.steps =  self.planetPositions[0,0,:].size




        self.planetPosFunction = inter.interp1d(self.time, self.planetPositions)


        self.calc_trad_parameters(self.destination_planet)
        self.angle_between_planets -=0.02
        self.calc_time_to_burn(self.angle_between_planets)

        self.launch_angle = self.calc_tangental_angle(self.time_to_launch)


        if (self.use_time0):
            self.time_to_launch = 0

        print "e home: ", self.system.e[0]
        print "e destinaton: ", self.system.e[self.destination_planet]



        print "Angle: ",self.calc_tangental_angle(self.time_to_launch)

        self.optimal_dist = self.calc_optimal_dist()
        print "Optimal radius of orbit: ",self.optimal_dist


        self.position = np.array([self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0]]) +np.array([-np.sin(self.calc_tangental_angle(self.time_to_launch)),np.cos(self.calc_tangental_angle(self.time_to_launch))])*(self.system.radius[self.home_planet]*self.km_to_au)

        v0 = [-np.sin(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .894*self.calc_influence()),np.cos(self.calc_tangental_angle(self.time_to_launch))*(1.0*self.dv_mainburn + .894*self.calc_influence())]





        self.v0_save = v0

        if not self.use_time0:
            self.planet_vel = self.calc_vel(0,self.time_to_launch)
        else:
            self.planet_vel = np.array([self.system.vx0[0],self.system.vy0[0]])

        print norm(self.planetPosFunction(self.time_to_launch)[:,0] - self.position)/self.km_to_au
        print "Starting Pos: ", self.position
        print "Total speed: ",norm(np.array(v0))





        self.velocity = np.array(v0) + self.planet_vel
        self.pos_over_time[:,0] = self.position





    def main_sequence(self):

        self.e.startEngine()
        self.force_per_box = self.e.force_gained
        self.mass_lost_per_box_per_sec = self.e.mass_escaped_per_sec
        self.force_engine = self.force_per_box*self.numb_boxes
        self.mass_lost_per_sec = self.mass_lost_per_box_per_sec*self.numb_boxes
        print "Engine ignited"
        self.fuel_cal(self.delta_v)
        print "Lauching with ", self.fuel, " kg fuel."
        print "---------------------"


    def fuel_cal(self,delta_v):
        self.fuel = self.e.launcher_mass*(np.exp(self.e.mass_escaped_per_sec*delta_v/self.e.force_gained) - 1)


    def solarcellSize(self,power):
        eff = 0.12
        return 1.0/eff * power*((self.starRadius*self.km_to_au)**2/(self.system.a[1])**2)*self.starTemp

    def acceleration(self,planetPos):


        r = planetPos[:,:] - self.position[:,np.newaxis]

        return np.sum(self.G*self.planetMasses[:]/(norm(r,axis = 0)**3)*r,axis = 1) - (self.G*self.starMass/(norm(self.position,axis = 0)**3))*self.position

    def main_loop(self):


        print "---------------------------------"
        print "Beginning to calculate the orbit: "

        #time = np.zeros(self.steps_per_year*self.sim_time)
        if (self.use_time0):
            time = 0.
        else:
            time =self.time_to_launch

        time_save = [time]

        #print self.planetPosFunction(time[-1])

        r_dest = (self.calc_dist_to_planet(self.destination_planet,time,self.position))
        r_home = (self.calc_dist_to_planet(self.home_planet,time,self.position))

        min_dist_to_planet = 200
        time_closest_incounter = 0
        index_closest_incounter = 0

        dt_mod = 1000
        dt = self.dt/dt_mod

        start = tid.clock()

        times_of_burns = [time]

        couter = 0
        r_dest = 100

        lauched = False
        close_to_dest = False
        inject_burn = False
        correct_burn = False
        check_feq = 1000

        at_correct_data = np.zeros(4)
        at_inject_data = np.zeros(4)

        check1 = False
        check2 = False
        check1_time = time + 1
        check2_time = time + 2

        check1_data = np.zeros(5) #pos,vel,time
        check2_data = np.zeros(5)



        self.velocity -= 0.5*self.acceleration(self.planetPosFunction(time))*self.dt/dt_mod

        for i in xrange(1,int(self.steps_per_year*self.sim_time)):


            if (i%check_feq == 0):
                r_dest = (self.calc_dist_to_planet(self.destination_planet,time,self.position))
                r_home = (self.calc_dist_to_planet(self.home_planet,time,self.position))

                if r_home > 0.005 and not lauched:

                    print "You have lauched"
                    lauched = True
                elif r_dest < 0.0005 and not close_to_dest:
                    print "Getting close to destinaton"
                    dt_mod = 1000
                    close_to_dest = True
                elif r_dest > 0.0005 and close_to_dest:
                    print "Getting away from the planet"
                    close_to_dest = False
                elif r_dest < 0.003 and not correct_burn:
                    at_correct_data[:2] = self.position
                    at_correct_data[2:] = self.velocity
                    print "Doing a corretion burn at time ",time
                    self.correction_burn = self.calc_correction_burn(time,.5)#.61)
                    print "With velocity ", self.correction_burn
                    self.velocity += self.correction_burn
                    correct_burn = True
                    check_feq = 10
                    times_of_burns.append(time)
                elif r_dest < self.optimal_dist and not inject_burn:
                    at_inject_data[:2] = self.position
                    at_inject_data[2:] = self.velocity
                    self.injection_burn = self.calc_injection_burn(time)
                    self.velocity += self.injection_burn
                    print "Burning for orbit, at time ", time
                    print "With velocity ", self.injection_burn
                    inject_burn = True
                    times_of_burns.append(time)


                if time > check1_time and not check1:
                    check1_data[0] = self.position[0]
                    check1_data[1] = self.position[1]
                    check1_data[2] = self.velocity[0]
                    check1_data[3] = self.velocity[1]
                    check1_data[4] = time
                    check1 = True
                elif time > check2_time and not check2:
                    check2_data[0] = self.position[0]
                    check2_data[1] = self.position[1]
                    check2_data[2] = self.velocity[0]
                    check2_data[3] = self.velocity[1]
                    check2_data[4] = time
                    check2 = True


            if (close_to_dest or not lauched):
                dt = self.dt/dt_mod
                for j in xrange(dt_mod):
                    self.velocity += self.acceleration(self.planetPosFunction(time + j*dt))*dt
                    self.position += self.velocity*dt

            else:
                dt = self.dt
                self.velocity += self.acceleration(self.planetPosFunction(time + dt))*dt
                self.position += self.velocity*dt


            time += self.dt





            if (r_dest < min_dist_to_planet):
                min_dist_to_planet = r_dest
                time_closest_incounter = time
                index_closest_incounter = float(i)/self.writeingFreq




            if ((i)%self.writeingFreq == 0 ):
                self.pos_over_time[:,i/self.writeingFreq] = self.position
                time_save.append(time)









        print ""


        print "It took", (tid.clock()-start), " sec"
        #time_save = time



        print "Shortest distance to planet: ", min_dist_to_planet
        print "At time: ",time_closest_incounter
        print "Planet position is: ", self.planetPosFunction(time_closest_incounter)[:,self.destination_planet]
        print "Satelite position is: ",self.pos_over_time[:,int(index_closest_incounter)]


        if (inject_burn and correct_burn):
            print "-----------------------"
            print "dv and Fuel"
            print "Hohmann burn: ", norm(self.dv_mainburn)
            print "Main Burn: ", norm(self.v0_save)
            print "Correction Burn: ",norm(self.correction_burn)
            print "Injection Burn: ", norm(self.injection_burn)
            print "Total: ", norm(self.v0_save) + norm(self.correction_burn) +norm(self.injection_burn)
            print "Total fuel: ", self.calc_fuel(norm(self.v0_save)+norm(self.correction_burn)+norm(self.injection_burn))
            print "Fuel with extra: ", self.calc_fuel(norm(self.v0_save)+norm(self.correction_burn)+norm(self.injection_burn))*1.1

        self.dump_to_file(inject_burn,correct_burn,times_of_burns,check1_data,check2_data,at_correct_data,at_inject_data)


        plt.plot(self.pos_over_time[0,0],self.pos_over_time[1,0],"x")
        plt.plot(self.pos_over_time[0,:],self.pos_over_time[1,:])
        plt.plot(self.pos_over_time[0,-1],self.pos_over_time[1,-1],"o")

        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0],"rv")
        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0],self.planetPosFunction(self.time_to_launch)[1,0],"gv")
        plt.plot(self.planetPosFunction(self.time_to_launch)[0,0] - np.sin(self.calc_tangental_angle(self.time_to_launch))*0.01,self.planetPosFunction(self.time_to_launch)[1,0] + np.cos(self.calc_tangental_angle(self.time_to_launch))*0.01,"yv")

        pos_planets =  self.planetPosFunction(np.array(time_save))

        for p in xrange(self.numberOfPlanets):
            plt.plot(pos_planets[0,p,:],pos_planets[1,p,:])


        plt.plot(self.planetPosFunction(time_closest_incounter)[0,self.destination_planet],self.planetPosFunction(time_closest_incounter)[1,self.destination_planet],"*")



        r_relative_from_dest = self.planetPosFunction(np.array(time_save))[:,self.destination_planet] - self.pos_over_time




        plt.axis("equal")
        plt.xlabel("Position [AU]")
        plt.ylabel("Position [AU]")
        plt.title("Jourey of the Satellite")
        plt.show()

        plt.plot(r_relative_from_dest[0],r_relative_from_dest[1])
        plt.plot(0,0,"ro")
        plt.axis("equal")
        plt.xlabel("Position [AU]")
        plt.ylabel("Position [AU]")
        plt.title("The Satellite Relative to Isskji")
        plt.show()

        if self.save:
            np.save("posOverTime.npy",self.pos_over_time)
            np.save("time.npy",time_save)


    def calc_escape(self):

        return np.sqrt(2*self.G*self.system.mass[self.home_planet]/(self.system.radius[self.home_planet]*self.km_to_au))

    def calc_influence(self):
        k = 10
        r = norm(np.array([self.system.x0[0],self.system.y0[0]]))
        r_soi = r*(self.planetMasses[0]/self.starMass)**(2./5)

        r_inf = r/np.sqrt(k*self.starMass/self.planetMasses[0])
        return np.sqrt(-(2*self.G*self.planetMasses[0])/(r_soi) + 2*self.G*self.planetMasses[0]/(self.system.radius[0]*self.km_to_au))

    def calc_trad_parameters(self,planet_number):

        time_mod = -0.

        r1 = norm(np.array([self.system.x0[0],self.system.y0[0]]))
        r2 = norm(np.array([self.system.x0[planet_number],self.system.y0[planet_number]]))


        self.v_mainburn = np.sqrt(2*(self.G*self.starMass/r1 - self.G*self.starMass/(r1+r2)))
        self.dv_mainburn = np.sqrt(self.G*self.starMass/r1)*(np.sqrt(2*r2/(r1+r2)) - 1)
        self.time_to_encounter = time_mod + np.pi*np.sqrt((r1+r2)**3/(8*self.G*self.starMass))

        self.omega_destination = 2*np.pi/(np.sqrt((4*np.pi**2 * self.system.a[planet_number]**3)/(self.G*(self.starMass+self.planetMasses[planet_number]))))
        self.omega_home = 2*np.pi/(np.sqrt((4*np.pi**2 * self.system.a[self.home_planet]**3)/(self.G*(self.starMass+self.planetMasses[self.home_planet]))))
        self.angle_between_planets = np.pi - self.omega_destination*self.time_to_encounter

        print "Burn required: ",self.dv_mainburn
        print "Trip takes: ", self.time_to_encounter
        print "Planet moves at angular v: ", self.omega_destination
        print "Home planets anguar v: ", self.omega_home
        print "Angle between planet when burn: ",self.angle_between_planets
        print "Escape Velocity for planet: ",self.calc_escape()

    def calc_time_to_burn(self,angle):

        print "Beginning to calculate burn time. This may take a couple of minutes..."

        min_angle = 10000.
        time_for_min_angle = 0
        #time = 0.0
        time = 6.
        eps = 1.e-4

        sim_time = 20.
        steps = 20000
        dt = 1./steps
        for i in xrange(1,int(steps*sim_time)-1):
            angle_home = np.arctan2(self.planetPosFunction(time)[1,self.home_planet],self.planetPosFunction(time)[0,self.home_planet])
            angle_destination = np.arctan2(self.planetPosFunction(time)[1,self.destination_planet],self.planetPosFunction(time)[0,self.destination_planet])

            time += dt

            if abs((angle_destination - angle_home)-angle) < eps:
                break

            print (float(i)/(steps*sim_time))*100, "%            \r",
        print ""


        self.time_to_launch = time
        print "Smallest difference in angle: ",abs((angle_destination - angle_home)-angle)
        print "At time: ",self.time_to_launch


    def calc_dist_to_planet(self,planet_number,time,pos):

        r = norm(self.planetPosFunction(time)[:,planet_number] - pos)
        return abs(r)

    def calc_tangental_angle(self,time):
        return  np.arctan2(self.planetPosFunction(time)[1,self.home_planet],self.planetPosFunction(time)[0,self.home_planet])


    def calc_vel(self,planet_number,time):
        dx = (self.planetPosFunction(time+self.dt)[0,planet_number] - self.planetPosFunction(time-self.dt)[0,planet_number])/(2*self.dt)
        dy = (self.planetPosFunction(time+self.dt)[1,planet_number] - self.planetPosFunction(time-self.dt)[1,planet_number])/(2*self.dt)

        return np.array([dx,dy])

    def calc_injection_burn(self,time):
        r =   self.position - self.planetPosFunction(time)[:,self.destination_planet]

        theta = np.arctan2(r[1],r[0])

        orb_vel = np.sqrt(self.G*self.planetMasses[self.destination_planet]/norm(r))

        return np.array([-orb_vel*np.sin(theta),orb_vel*np.cos(theta)])-self.velocity+self.calc_vel(self.destination_planet,time)

    def calc_fuel(self,dv):
        n_e = 2.120316e-13
        f_b = 1.72650965887e-09

        return 1100*(np.exp((dv*4743.7173611)*n_e/f_b)-1)

    def calc_optimal_dist(self):
        return self.system.a[self.destination_planet]*np.sqrt(self.planetMasses[self.destination_planet]/(10*self.starMass))

    def calc_correction_burn(self,time,factor):
        r = self.position - self.planetPosFunction(time+0.001)[:,self.destination_planet]
        normal_r = 1.0/(norm(r)) * r
        print "Length normal r: ", norm(normal_r)
        return -factor*normal_r

    def dump_to_file(self,inject,correct,times_of_burns,check1_data,check2_data,at_correct_data,at_inject_data):
        dv_from_orbit = 0
        outFile = open("data.txt","w")
        outFile.write("Starting Position: " + str(self.pos_over_time[:,0]) + "\n")
        outFile.write("V0: " + "Vel: " + str(self.v0_save) + "; Speed: " + str(norm(self.v0_save)) + " At time: " + str(times_of_burns[0]) + "\n")
        if (inject and correct):
            outFile.write("When correcting: Pos: " + str(at_correct_data[:2]) + "; Vel: " + str(at_correct_data[2:]) + "\n")
            outFile.write("Correction: " + "Vel: " + str(self.correction_burn) + "; Speed: " + str(norm(self.correction_burn)) + " At time: " + str(times_of_burns[1]) + "\n")
            outFile.write("When injecting: Pos: " + str(at_inject_data[:2]) + "; Vel: " + str(at_inject_data[2:]) + "\n")
            outFile.write("Injection: " + "Vel: " + str(self.injection_burn) + "; Speed: " + str(norm(self.injection_burn)) + " At time: " + str(times_of_burns[2]) + "\n")
            dv_from_orbit = norm(self.correction_burn) + norm(self.injection_burn)
        outFile.write("Total dv: " + str(norm(self.v0_save) + dv_from_orbit) + "\n")
        outFile.write("Fuel: " + str(self.calc_fuel(dv_from_orbit + norm(self.v0_save))) + "\n")
        outFile.write("Fuel with extra: " + str( self.calc_fuel( (dv_from_orbit + norm(self.v0_save))*1.1 ) ) + "\n")
        outFile.write("-----------------------------" + "\n")
        outFile.write("Check 1: pos " + str(check1_data[:2]) + "; Vel: " + str(check1_data[2:4]) + "; At time: " + str(check1_data[-1]) + "\n")
        outFile.write("Check 2: pos " + str(check2_data[:2]) + "; Vel: " + str(check2_data[2:4]) + "; At time: " + str(check2_data[-1]) + "\n")
        outFile.close()

destinaton_planet = 1
sat = satelite(destinaton_planet)
print "The size of the solar panels: ", sat.solarcellSize(40)
sat.main_loop()
