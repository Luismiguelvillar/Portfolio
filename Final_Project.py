"""

Theory Computing Project: Monte Carlo Simulation of Nuclear Reactor
-------------------------------------------------------------------

This code will produce a random walk of a number N_neutrons of neutrons through a 
spherical material. This material consists on a mixture of two kinds of Uranium
and a moderator, in this case H2O. The goal of this simulation is to calculate the 
multiplication factor of the reactor given certain know parameters, such as the 
composition of the reactor and its radius.

Code by
Luis Miguel Villar Padruno.

"""

import matplotlib as mpl
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.constants import Avogadro

""" Choose if you want or not to optimize """

low_optimization = True # Low Optimization removes the dead neutrons from the neutron list (DEFAULT = TRUE)
large_optimization = True # Large Optimization is used to reduce the number of neutrons when k is larger than 1. (DEFAULT = FALSE) (YA NO FUNCIONA)

""" Are you interested in delayed neutrons? """
delayed_analysis = True

#Are you interested in text for random walk and results of every iteration
text = True

# For thermal neutrons put 0, for fast neutrons put 1:
type_neutron = 1

# Define number of trials for generating neutrons in a sphere of radius R in cm.
N = 100000
R = 10 #cm

# Resolution is the increment taken into account when we determine if the sample is uniform.
resolution = 0.05

# This is the number of steps of the randomwalk.
steps = 10

# Default parts of U235, U238 and H20.
pH20 = 0.7
pU235 = 0.00216
pU238 = 0.29784

"""pH20 = 0
pU235 = 1
pU238 = 0"""

def validation_number(entry):
    """
    This function will validate that the entries are valid  
    """
    try:
        float(entry)
        return True
    except:
        return False

class Neutron:
    def __init__(self, position, name, father, n_fathers, dead, delayed, type):
        """
        Defines a class called Neutron that will contain the characteristic of each neutron of the simulation.
        """
        self.name = name
        self.position = position
        self.father = father
        self.n_fathers = n_fathers
        self.dead = dead
        self.delayed = delayed
        self.type = type

    def __repr__(self):
        """
        This function just tells the code how to represent the instances of the class Neutron.
        """
        return f"Neutron {self.name} at position ({self.position[-1][0]},{self.position[-1][1]},{self.position[-1][2]}) with father {self.father} and {self.n_fathers} fathers. DEAD = {self.dead}"

    def radio(self):
        """
        This function returns the radius at which the neutron currently is.
        """
        return np.sqrt(self.position[-1][0]**2 + self.position[-1][1]**2 + self.position[-1][2]**2)

def k_value(N_alive, N_alive1):
    """
    Calculates de k value using the total number of neutrons alive and the total
    number of neutrons alive a generation after.

    Parameters:
    -----------
    N_alive: Number of neutrons alive in the current generation.
    N_alive1: number of neutrons alive in the following generation.

    Return:
    -------
    k_value: float
    """
    return N_alive1 / N_alive

def prom_k(list_k):
    """
    This function receives a list of K's for each generation and calculates
    the average k over a certain number of generations. If the number of neutrons
    alive is zero for the next generation it does not take it into consideration.

    Parameters:
    -----------
    list_k: list of float k values.

    Return:
    -------
    float.
    """
    cont = 0
    suma = 0
    n_num = 0 
    for element in list_k:
        if element != "N" and cont > 1:
            suma = suma + element
            n_num = n_num + 1
        cont = cont + 1
    return suma / n_num

def parameters_of_randomwalk(parts_H20, parts_U235, parts_U238):
    """
    This functionm calculates the average mean free path of each neutron in the mixture and the probs
    of scattering, absorption and fission. 

    Parameters:
    -----------
    parts_H20: float representing the fraction of water in mixture.
    parts_U235: float representing the fraction of Uranium 235 in mixture.
    parts_U238: float representing the fraction of Uranium 238 in mixture.

    Return:
    -------
    parameters: list containing the average mean free path and the probabilities of absorption,
    sccattering and fission. The first entry of the list will be a list with the thermal neutron
    parameters and the second will be a list with the fast neutron parameters.

    """
    N_A = Avogadro
    rho_U235 = 18.7 # g/cm^3
    rho_U238 = 18.9 # g/cm^3
    rho_H20 = 1 # g/cm^3
    molar_mass_U235 = 235 #g/mol
    molar_mass_U238 = 238 #g/mol
    molar_mass_H20 = 18 #g/mol

    # Now we calculate the density of atoms for U235, U238 and H20
    N_U235 = ((rho_U235 * N_A) / molar_mass_U235) * parts_U235
    N_U238 = ((rho_U238 * N_A) / molar_mass_U238) * parts_U238
    N_H20 = ((rho_H20 * N_A) / molar_mass_H20) * parts_H20

    # Now we define the microscopical cross sections of each component for thermal neutrons:
    U235_abs_therm = 99 # barns
    U235_scat_therm = 10 # barns
    U235_fis_therm = 583 # barns
    U238_abs_therm = 2 # barns
    U238_scat_therm = 9 # barns
    U238_fis_therm = 0.00002 # barns
    H20_abs_therm = 0.66 # barns
    H20_scat_therm = 4 # barns

    # Now we do the same for fast neutrons:
    U235_abs_fast = 0.09 # barns
    U235_scat_fast = 4 # barns
    U235_fis_fast = 1 # barns
    U238_abs_fast = 0.07 # barns
    U238_scat_fast = 5 # barns
    U238_fis_fast = 0.3 # barns
    H20_abs_fast = 0.26 # barns
    H20_scat_fast = 4.44 # barns

    # Now we calculate the cross sections for the mixture (barns)
    abs_micro_sigma_thermal = parts_H20 * H20_abs_therm \
        + parts_U235 * U235_abs_therm + parts_U238 * U238_abs_therm
    abs_micro_sigma_fast = parts_H20 * H20_abs_fast \
          + parts_U235 * U235_abs_fast + parts_U238 * U238_abs_fast
    fis_micro_sigma_thermal = parts_U235 * U235_fis_therm \
          + parts_U238 * U238_fis_therm
    fis_micro_sigma_fast = parts_U235 * U235_fis_fast \
        + parts_U238 * U238_fis_fast
    scat_micro_sigma_thermal = parts_H20 * H20_scat_therm \
          + parts_U235 * U235_scat_therm + parts_U238 * U238_scat_therm
    scat_micro_sigma_fast = parts_H20 * H20_scat_fast \
        + parts_U235 * U235_scat_fast + parts_U238 * U238_scat_fast

    total_micro_thermal = abs_micro_sigma_thermal + \
        fis_micro_sigma_thermal + scat_micro_sigma_thermal
    total_micro_fast = abs_micro_sigma_fast + \
        fis_micro_sigma_fast + scat_micro_sigma_fast

    # We calculate the macroscopical cross sections, but the result will be in barn/cm^3. 
    # Therefore we multiply the result by 10^(-24) and the new unit is cm^-1
    macroscopical_sigma_thermal = (N_H20 * (H20_abs_therm + H20_scat_therm) + \
                                   N_U235 * (U235_abs_therm + U235_fis_therm + U235_scat_therm) + \
                                    N_U238 * (U238_abs_therm + U238_fis_therm + U238_scat_therm)) * 10**(-24)
    macroscopical_sigma_fast = (N_H20 * (H20_abs_fast +  H20_scat_fast) + \
                                N_U235 * (U235_abs_fast +  U235_fis_fast +  U235_scat_fast) + \
                                    N_U238 * (U238_abs_fast + U238_fis_fast + U238_scat_fast)) * 10**(-24)

    # Now we define the average mean free path as:
    lambda_thermal = 1 / macroscopical_sigma_thermal # cm
    lambda_fast = 1 / macroscopical_sigma_fast # cm

    # The probabilities of absorption, scatter and fission for thermal and fast neutrons are:
    probability_abs_thermal = abs_micro_sigma_thermal / total_micro_thermal
    probability_abs_fast = abs_micro_sigma_fast / total_micro_fast
    probability_scat_thermal = scat_micro_sigma_thermal / total_micro_thermal
    probability_scat_fast = scat_micro_sigma_fast / total_micro_fast
    probability_fis_thermal = fis_micro_sigma_thermal / total_micro_thermal
    probability_fis_fast = fis_micro_sigma_fast / total_micro_fast

    return [[lambda_thermal, probability_abs_thermal, probability_fis_thermal, probability_scat_thermal ],
            [lambda_fast, probability_abs_fast, probability_fis_fast, probability_scat_fast]]

def neutrones(sph, num_neutrons):
    """
    The function takes a list of positions and the number of initial neutrons and
    returns a list of Neutron instances (objects).

    Parameters:
    ----------
    sph: list of list representing positions of neutrons.
    num_neutrons: int representing the initial number of neutrons.

    Return:
    -------
    neutro: list of Neutron instances (objects).

    """
    neutro = []

    for i in range(0,num_neutrons):
        temp = Neutron([list(np.delete(sph[i],-1, axis=0))], str(i), None, 0, 0, False, type_neutron)
        neutro.append(temp) 
    return neutro

def prob(prob):
    """
    This functions takes a probability and then generates a random number from 0 to 1.
    Then it returns True if the number is less or equal to the probability we want, and
    False otherwise.

    Parameters:
    -----------
    prob: float from 0 to 1, representing the probability of success.

    Return:
    -------
    Bool: False of True of event occurring.

    """
    number = rand.uniform(low=0, high = 1)
    if number <= prob and number != 0:
        return True
    else: 
        return False

def sample_cube(number_of_trials, R):
    """
    This function receives an integer and returns a uniform sample of points in a square.
    The return object is a list of lists, which consists of three entries, each is a list
    of the copordinates of a point in the space.

    Parameters:
    -----------
    number_of_trials: number of points generated in cube.
    R: float representing the radius of the sphere desired (side of the cube over 2).

    Return:
    -------
    cube: list of positions in the form list of list of 3 coordinates.

    """

    cube = rand.uniform(low = -R, high = R, size = (number_of_trials, 3))
    x_cube = cube[:,0]
    y_cube = cube[:,1]
    z_cube = cube[:,2]
    return cube

def sample_sphere(cube, R):
    """
    This function takes a list of lists containming the positions of neutrons
    in the cube of side R and returns a list with the indices of the elements
    in cube that meet the condition x^2 + y^2 + z^2 <= R.

    Parameters:
    -----------
    cube: list of the form [Position_0, Position_1, ...] & Position_i => [x_i,y_i,z_i]
    R: float representing the radius of the desired sphere.

    Return:
    -------
    in_sphere: list of positions of neutrons that fulfill the condition above.

    """
    x_cube = cube[:,0]
    y_cube = cube[:,1]
    z_cube = cube[:,2]
    cube = [x_cube, y_cube, z_cube]

    in_sphere = np.where(np.sqrt(cube[0]**2 + cube[1]**2 + cube[2]**2) <= R)
    return in_sphere

def directions(N_neutrons):
    """
    This function will create random directions for each neutron generated within the
    sphere. It receives an integer and will return a list of uniform directions in cartesian coordinates. 
    All the vectors are normalized.

    Parameters:
    -----------
    N_neutrons: int number of neutrons in sphere.

    Return:
    -------
    directions: list of list of direction vectors of the form [x, y, z]

    """
    directions = rand.uniform(low=-1, high=1, size=(N_neutrons, 3))
    directions_normalized = []
    factors = []
    for element in directions:
        factors.append(np.sqrt(element[0]**2 + element[1]**2 + element[2]**2))
    contador = 0
    for element in directions:
        directions_normalized.append(element / factors[contador])
        contador = contador + 1

    return np.array(directions_normalized)

def graph(cube, sphere, N_neutrons):
    """
    This function will receive two lists and then produce a scatter plot of the neutrons in
    a sphere of radius R.

    Parameters:
    -----------
    cube: list of points in a cube of side 2R
    sphere: list of indeces of the points that are inside the sphere.

    Return:
    -------
    None: (But shows a graph as an emergent window).

    """
    x_cube = cube[:,0]
    y_cube = cube[:,1]
    z_cube = cube[:,2]
    cube = [x_cube, y_cube, z_cube]

    x_in_sphere = cube[0][sphere[0]]
    y_in_sphere = cube[1][sphere[0]]
    z_in_sphere = cube[2][sphere[0]]

    mpl.rcParams["legend.fontsize"] = 10

    fig = plt.figure()
    ax = plt.axes(projection = "3d")
    ax.scatter(x_in_sphere, y_in_sphere, z_in_sphere, s = 1, color = "green")
    ax.set_aspect('equal')
    ax.set_ylabel("y axis [cm]")
    ax.set_xlabel("x axis [cm]")
    ax.set_zlabel("z axis [cm]")
    ax.grid("on")

    plt.title(label = str(N_neutrons) + " Neutrons in Sphere of R = " + str(R) + "cm", fontsize = 20 )
    plt.savefig("NeutronsInSphere.png", dpi = 1000)
    plt.show()

def sphere_points(cube, sphere):
    """
    This fuinction receives the cube random sample and the sphere indices and then returns
    the points of the sphere in the format [Position_1, Position_2, ...].

    Parameters:
    -----------
    cube: list of the cube random sample in the usual form [Position_1, Position_2, ...]. 
    sphere: list containing the indices of the points that are inside a sphere of radius R.

    Return:
    -------
    sphere: list of list of the form [Position_1, Position_2, ...] with the positions of the
    neutrons in the sphere of radius R.

    """
    sphere = cube[sphere[0]]
    new_sphere = []
    cont = 0
    for neutron in sphere:
        neutron = list(neutron)
        neutron = np.append(neutron, int(cont))
        new_sphere.append(neutron)
        cont = cont + 1
    sphere = new_sphere

    return sphere

def isUniform(sphere_points, resolution, ranwalke):
    """
    This function receives a list of points and determines if the neutrons are uniformily distributed in the sphere.
    The function will determine the distance of each neutron to the centre of the sphere and then it will create a list
    called cumulative which will have the number of neutrons presented at a radius < R all of this scaled by 1/r^3.

    Parameters:
    -----------
    sphere_points: list of points in the usual form [Position_1, Position_2, ...].
    resolution: float representing the value of the increment when performing the calculation.

    Return:
    -------
    None (Graph showing relation).

    """
    if ranwalke == None:
        distances_centre = []

        for i in range(len(sphere_points)):
            distance = np.sqrt(sphere_points[i][0]**2 + sphere_points[i][1]**2 + sphere_points[i][2]**2)
            distances_centre.append(distance)
        
        comulative = []
        radius_list = []
        radius = resolution
        while radius < R:
            N_within = 0
            for element in distances_centre:
                if element < radius:
                    N_within = N_within + 1
            radius_list.append(radius)
            comulative.append(N_within/radius**3)
            radius = radius + resolution

        fig = plt.figure()
        ax = plt.subplot(111)

        ax.scatter(radius_list, comulative, s = 0.8, color = "red")
        ax.set_xlabel("Radius" + r'$[cm]$')
        ax.set_ylabel("Total Number of neutrons / " + r'$r^{3}$')
        plt.title("Total number of neutrons vs 1 / " + r'$r^{3}$')
        plt.grid('on')
        plt.savefig("Uniformity1.png", dpi = 1000)
        plt.show()

    else:
        final_pos = []
        for neutron in ranwalke:
            if neutron.dead == False:
                final_pos.append(neutron.position[-1])

        distances_centre = []

        for i in range(len(final_pos)):
            distance = np.sqrt(final_pos[i][0]**2 + final_pos[i][1]**2 + final_pos[i][2]**2)
            distances_centre.append(distance)
        
        comulative = []
        radius_list = []
        radius = resolution
        while radius < R:
            N_within = 0
            for element in distances_centre:
                if element < radius:
                    N_within = N_within + 1
            radius_list.append(radius)
            comulative.append(N_within/radius**3)
            radius = radius + resolution

        fig = plt.figure()
        ax = plt.subplot(111)

        ax.scatter(radius_list, comulative, s = 0.8, color = "blue")
        ax.set_xlabel("Radius" + r'$[cm]$')
        ax.set_ylabel("Total Number of neutrons / " + r'$r^{3}$')
        plt.title("Total number of neutrons vs 1 / " + r'$r^{3}$' + "after " + str(steps) +" steps")
        plt.grid('on')
        plt.savefig("UniformityAfter.png", dpi = 1000)
        plt.show()

    return None

def N_delayed(neutrons):
    """
    This function receives a list of neutrons and then determine the number of delayed neutrons and then
    counts how many neutrons are alive and delayed.

    Parameters:
    -----------
    neutrons: list of instances of the class Neutron.

    Return:
    -------
    alive_delayed: Number of delayed neutrons (int).
 
    """
    alive_delayed = 0
    delayed = 0
    for n in neutrons:
        if n.dead == False and n.delayed == True:
            alive_delayed = alive_delayed + 1
        if n.delayed == True:
            delayed = delayed + 1

    return alive_delayed

def randomWalk(neutrons, steps, parameters):
    """
    Try of random Walk in 3D for N particles. This function receives a list of Neutron instances and the
    desired number of steps. Then performs a random walk of the neutrons, in each step the neutrons can
    be either aborbed, scattered or produce a Uranium atom fission, creating 2 new 'child' neutrons.

    Parameters:
    -----------
    neutrons: list of Neutron instances of the form [Neutron_1, Neutron_2, ...].
    steps: int value representing the desired number of steps for the random walk.

    Return:
    -------
    ranwalk: list of neutros again of the form [Neutron_1, Neutron_2, ...] but with 
    different attributes for the instances.

    """
    total_fisiones = 0
    total_abs = 0
    total_escape = 0
    k_calculated = []

    if text == True:
        print("The step 0 has: "+ str(len(neutrons)) + " Neutrons.")

    for step in range(0, steps):
        if text == True:
            print("This is the step:", str(step + 1), "out of", str(steps))
            print("---")
        init_number = len(neutrons)

        fis_n = 0
        abs_n = 0
        esc_n = 0
        alive_n = 0 
        alive_n1 = 0
        delayed_produced = 0
        fissions_by_delayed = 0

        for n in neutrons:
            if n.dead == 0:
                alive_n = alive_n + 1

        cont = 0
        for neutron in neutrons:
            temp = directions(1)[0] * parameters[neutron.type][0]
            new_position = [temp[0] + neutron.position[-1][0], temp[1] + neutron.position[-1][1] , temp[2] + neutron.position[-1][2]]

            # Validando fuera del reactor
            if np.sqrt(new_position[0]**2 + new_position[1]**2 + new_position[2]**2) < R and neutron.dead == 0 and cont < init_number:
 
                # Validando si ocurre fision
                if prob(parameters[neutron.type][2]) == True:
                    if neutron.delayed == True:
                        fissions_by_delayed = fissions_by_delayed + 1

                    temp1 = directions(1)[0] * parameters[1][0]
                    temp2 = directions(1)[0] * parameters[1][0]
                    new_position1 = [temp1[0] + neutron.position[-1][0], temp1[1] + neutron.position[-1][1] , temp1[2] + neutron.position[-1][2]]
                    new_position2 = [temp2[0] + neutron.position[-1][0], temp2[1] + neutron.position[-1][1] , temp2[2] + neutron.position[-1][2]]
                    if neutron.type == 1:
                        neutrons.append(Neutron([new_position1], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))
                        neutrons.append(Neutron([new_position2], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))
                    if neutron.type == 0:
                        neutrons.append(Neutron([new_position1], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))
                        neutrons.append(Neutron([new_position2], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))
                    neutron.dead = 1
                    if prob(0.42) == True:
                        temp3 = directions(1)[0] * parameters[1][0]
                        new_position3 = [temp3[0] + neutron.position[-1][0], temp3[1] + neutron.position[-1][1] , temp3[2] + neutron.position[-1][2]]
                        if neutron.type == 1:
                            neutrons.append(Neutron([new_position3], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))
                        if neutron.type == 0:
                            neutrons.append(Neutron([new_position3], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, False, 1))

                    if prob(0.0065) == True and delayed_analysis == True:
                        temp4 = directions(1)[0] * parameters[neutron.type][0]
                        new_position4 = [temp4[0] + neutron.position[-1][0], temp4[1] + neutron.position[-1][1] , temp4[2] + neutron.position[-1][2]]
                        neutrons.append(Neutron([new_position4], str(len(neutrons)), neutron.name, neutron.n_fathers + 1, 0, True, 0))
                        delayed_produced = delayed_produced + 1

                    #print("The neutron", neutron.name, "provocated a fission in step:", str(step))
                    total_fisiones = total_fisiones + 1
                    fis_n = fis_n + 1

                # Validando si ocurre abs
                elif prob(parameters[neutron.type][1]) == True and neutron.dead == 0:
                    total_abs = total_abs + 1
                    abs_n = abs_n + 1
                    #print("The neutron", neutron.name, "was absorbed in step:", str(step))
                    neutron.dead = 1

                else:
                    neutron.position.append(new_position)

            elif np.sqrt(new_position[0]**2 + new_position[1]**2 + new_position[2]**2) >= R and neutron.dead == 0:
                #print("The neutron", neutron.name, "escaped in step:", str(step))
                #print("The radius is:", str(neutron.radio()) ,"and the new radius would be:", str(np.sqrt(new_position[0]**2 + new_position[1]**2 + new_position[2]**2)))
                total_escape = total_escape + 1
                esc_n = esc_n + 1
                neutron.dead = 1

            else:
                #print("NO CONTABILIZA")
                #print("the radius is:", str(neutron.radio()))
                #print("the counter is:",  str(cont), "and init_number is:", str(init_number))
                #print("Is dead:", neutron.dead)
                None
            cont = cont + 1

        for n in neutrons:
            if n.dead == 0:
                alive_n1 = alive_n1 + 1

        # Calculate K:
        if alive_n != 0 and alive_n1 !=0:
            k_step = k_value(alive_n, alive_n1)
        else:
            k_step = "N"

        k_calculated.append(k_step)
        

        # Contar delayed neutrons
        if delayed_analysis:
            delayed = N_delayed(neutrons)
            if text == True:
                print("Delayed neutrons alive:", delayed)
                print("Delayed neutrons produces in step:", delayed_produced)
                print("---")

        # OPTIMIZAR
        
        if low_optimization:
            neutron_temp = []
            for ntr in neutrons:
                if ntr.dead == 0:
                    neutron_temp.append(ntr)
            neutrons = neutron_temp

        #LARGE OPTIMIZATION
        if large_optimization:
            while len(neutrons) >  init_number:
                aux = 0
                for neutron in neutrons:
                    if len(neutrons) > init_number:
                        if prob(0.5) == True:
                            neutrons.pop(aux)
                    aux = aux + 1
        if text == True:
            print("K: " + str(k_step))

    N_alive = 0
    for neut in neutrons:
        if neut.dead == False:
            N_alive = N_alive + 1
    
    return neutrons, total_abs, total_fisiones, N_alive, total_escape, k_calculated, fissions_by_delayed

def graph_randonWalk(ranwalk):

    """
    This function takes a list ranwalk of neutrons and generates a plot describing the trayectory of each neutron
    through the reactor.

    Parameters:
    -----------
    ranwalk: list of neutron instances of the form [Neutron_1, Neutron_2, ...].

    Return:
    -------
    None (Graph showing paths of each Neutron)

    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for neutron in ranwalk:
        if neutron.dead == False:
            n_x = [pos[0] for pos in neutron.position]
            n_y = [pos[1] for pos in neutron.position]
            n_z = [pos[2] for pos in neutron.position]
            ax.plot(n_x, n_y, n_z)

    ax.set_xlabel("x axis [cm]")
    ax.set_ylabel("y axis [cm]")
    ax.set_zlabel("z axis [cm]")
    plt.title(label = "Trajectories of the living neutrons inside the reactor")
    plt.savefig("Trajectories.png", dpi = 1000)
    plt.show()

def graph_final(ranwalk, N_alive):
    """
    This function will receive the neutrons of the random walk and the produce a scattered
    graph of the final positions of the neutrons that are alive.

    Parameters:
    -----------
    ranwalk: list of neutrons

    Return:
    -------
    None (Graph)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    final_pos = []
    for neutron in ranwalk:
        if neutron.dead == False:
            final_pos.append(neutron.position[-1])
    
    nx = [pos[0] for pos in final_pos]
    ny = [pos[1] for pos in final_pos]
    nz = [pos[2] for pos in final_pos]

    ax.scatter(nx, ny, nz, s = 1)
    ax.set_aspect('equal')
    ax.set_xlabel("x axis [cm]")
    ax.set_ylabel("y axis [cm]")
    ax.set_zlabel("z axis [cm]")
    plt.title(label = "Final positions of the " + str(N_alive) + " remaining neutrons in the reactor", fontsize = 20 )
    plt.savefig("FinalNeutrons.png", dpi = 1000)
    plt.show()

def graph_k(k_values, average_k, standard_deviation):
    """
    This function takes a list of values of k and then it plots them against the number of generations.
    
    Parameters:
    -----------
    values_k: list of floats representing calculated values of the multiplication factor k.

    Return:
    -------
    None (Graph)
    """

    new = []
    for element in k_values:
        if element != "N":
            new.append(element)
        if element == "N":
            new.append(None)

    fig, ax = plt.subplots()

    cont = 0
    new2 = []
    for i in new:
        new2.append(cont)
        cont = cont + 1

    ax.plot(new)
    ax.scatter(new2, new)
    ax.axhline(y = average_k + standard_deviation, color = 'red', linestyle ='-', label = r'$+\sigma$')
    ax.axhline(y = average_k - standard_deviation, color = 'red', linestyle ='-', label = r'$-\sigma$')
    ax.set_xlabel("Step")
    ax.set_ylabel("Value of K")
    plt.grid('on')
    plt.title("Measured Values of k against step number")
    plt.legend()
    plt.savefig("K_Value.png", dpi = 1000)
    plt.show()

def std_k(values_k):

    """
    This function calculates the standard deviation of the values measured for k,
    excluding bthe first 5, due to stabilization.
    """
    final_values_k = []
    cont = 0
    for element in values_k:
        if element != "N" and cont > 5:
            final_values_k.append(element)
        cont = cont + 1
    
    return np.std(final_values_k)

def validation():
    validation = False
    while validation == False:
        pH20 = input("Parts of water: ")
        pU235 = input("Parts of Uranium-235: ")
        pU238 = input("Parts of Uranium-238: ")

        if validation_number(pH20) ==  True and validation_number(pU235) == True and validation_number(pU238) == True:
            if float(pU235) >= 0 and float(pH20) >= 0 and float(pU238) >= 0:
                if float(pU238) + float(pH20) + float(pU235) == 1:
                    pH20 = float(pH20)
                    pU235 = float(pU235)
                    pU238 = float(pU238)
                    validation = True
                else:
                    print("The parts given do not sum up to 1.")
            else:
                print("One or more number given are negative, try again with positive numbers.")
        else:
            print("Inputs are not numbers, try again.")

    validation = False
    while validation == False:
        steps = input("The number of steps for the simulation: ")

        if validation_number(steps):
            if int(steps) > 5:
                steps = int(steps)
                validation = True
            else: 
                print("The number of steps should be greater than 5")

    validation = False
    while validation == False:
        R = input("Radius of the sphere in centimeters: ")
        if validation_number(R) == True:
            if float(R) >= 0:
                if float(R) > 0:
                    R = float(R)
                    validation = True
                else:
                    print("Radius given should be greater than 0.")
            else:
                print("Radius given is negative, try with a positive number.")
        else:
            print("Radius given is not a number, try again.")
    return pH20, pU235, pU238, R, steps

def Principal():

    """ Variables """
    pH20, pU235, pU238, R, steps = validation()
    cube = sample_cube(N, R)
    sphere_indices = sample_sphere(cube, R)
    sphere = sphere_points(cube, sphere_indices)
    N_neutrons = len(sphere_indices[0])
    neutron_list = neutrones(sphere, N_neutrons)
    parameters = parameters_of_randomwalk(pH20, pU235, pU238)

    """ Determining Uniformity """
    isUniform(sphere, resolution, None)

    """ Random Walk """
    # These are the average Mean free path in centimeters and the probabilities of fission and absorption. 
    lamb = parameters[type_neutron][0]
    prob_fission = parameters[type_neutron][2]
    prob_abs = parameters[type_neutron][1]

    # These are the list of neutrons, the total number of fissions, number of neutrons alive and number of leaked neutrons.
    ranwalk, total_abs, total_fis, N_alive, total_escape, k_values, fis_delayed = randomWalk(neutron_list, steps, parameters)

    """ Calculating the average k """
    average_k = prom_k(k_values)
    standard_deviation = std_k(k_values)

    """ Presenting Results """
    if text == True:
        print("---------------------------------")
        print("The mean free path of the fast neutrons is:", str(parameters[1][0]), "centimeters")
        print("The mean free path of the thermal neutrons is:", str(parameters[0][0]), "centimeters")
        print("The probability of fission of fast neutrons is:", str(parameters[1][2]))
        print("The probability of fission of thermal neutrons is:", str(parameters[0][2]))
        print("The probability of absorption of fast neutrons is:", str(parameters[1][1]))
        print("The probability of absorption of thermal neutrons is:", str(parameters[0][1]))
        print("The probability of scattering of fast neutrons is:", str(1 - (parameters[1][2] + parameters[1][1])))
        print("The probability of scattering of thermal neutrons is:", str(1 - (parameters[0][2] + parameters[0][1])))
        print("The average value of k is:", str(average_k), "+/-", str(standard_deviation))
        print("The total number of fissions was:", str(total_fis))
        print("The total number of fissions provocated by delayed neutrons was:", str(fis_delayed))
        print("The total number of absorbtions was:", str(total_abs))
        print("The initial number of neutron was:", str(N_neutrons))
        print("The number of neutrons alive is:", str(N_alive))
        print("---------------------------------")

    # Graphs
    graph(cube,sphere_indices, N_neutrons)
    graph_randonWalk(ranwalk)
    graph_final(ranwalk, N_alive)
    graph_k(k_values, average_k, standard_deviation)
    isUniform(None, resolution, ranwalk)

    return average_k, standard_deviation

Principal()