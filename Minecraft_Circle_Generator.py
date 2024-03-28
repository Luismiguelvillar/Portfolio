"""
Circle Generator
----------------
This program creates a circle using the midpoint algorithm. (WORK IN PROGRESS)

by Luis Miguel Villar Padruno
2023
"""

import numpy as np
import matplotlib.pyplot as plt
import tkinter
import ttkbootstrap as ttk

"""
Debemos definir un plano coordenado donde el centro de la circunferencia es
siempre el punto (0,0)
"""
thin =  True # The circle should be thick or just one block thin?
d = 12  #diametro del circulo en numeros naturales

def distancia_centro(coor):
    """
    This function calculates the distance between a block and the origin of the
    reference frame (which coincides with the center of the circunference)

    input:
    ------
    coor: list of the form [coor(x), coor(y)]

    return:
    -------
    float which represents the module of the vector [x,y]

    """
    return np.sqrt(coor[0]**2 + coor[1]**2)

def plot_list(lista_bloques):
    """
    This function recibes a list of the form [[x1, y1], [x2, y2], ..., [xn, yn]] and then builds
    an scatter plot with the points of the list.

    Input:
    ------
    lista_bloques: list of lists containing the coordenates of points in the space.

    Return:
    -------
    None (Plot)

    """
    xs = []
    ys = []
    
    for i in lista_bloques:
        xs.append(i[0])
        ys.append(i[1])

    figure, axes = plt.subplots()
    plt.scatter(xs, ys, c='red')
    plt.xlim(int(-d/2) - int(d/4), int(d/2) + int(d/4))
    plt.ylim(int(-d/2) - int(d/4), int(d/2) + int(d/4))
    axes.set_aspect(1)

    ticks_frequency = 1
    x_ticks = np.arange(int(-d/2) - int(d/4), int(d/2) + int(d/4)+1, ticks_frequency)
    y_ticks = np.arange(int(-d/2) - int(d/4), int(d/2) + int(d/4)+1, ticks_frequency)
    axes.set_xticks(x_ticks[x_ticks])
    axes.set_yticks(y_ticks[y_ticks])

    axes.set_xticks(np.arange(int(-d/2) - int(d/4), int(-d/2) - int(d/4)+1), minor=True)
    axes.set_yticks(np.arange(int(-d/2) - int(d/4), int(-d/2) - int(d/4)+1), minor=True)

    axes.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

    plt.show() 

    return None

def radius_treatment():
    """
    This function just selects the appropiate value for the radius of the circunference given
    a diameter.

    Input:
    ------
    None (uses global variable d, which is an integer)

    Return:
    -------
    radio: integer.
    coordenada: list representing the position of the first block.

    """
    if d % 2 == 0:
        radio = d/2
        coordenada = [0.5, d/2 - 0.5]
    else:
        radio = d/2
        coordenada = [0, d/2 - 0.5]
    return radio, coordenada

def midpoint_algorithm(coordenada, radio):
    """
    This function builds 1/8 of a circle using the midpoint algorithm, the points that describe
    this arc of circunference are on the form [x, y] and the coordinates x and y always fulfill
    the following conditions: x <= y , x > 0 & y > 0.

    Input:
    ------
    coordenada: list which represents the position of the first block ([x0, y0]). 
    radio: integer representing the desired radius of the circunderence.

    Return:
    -------
    bloques: list of list in the form [[x1, y1], [x2, y2], ..., [xn, yn]] representing the positions
    of the blocks which describe the arc of circunference.

    """
    bloques = []
    if d%2 != 0:
        bloques.append([0, d/2 - 0.5])
    else:
        bloques.append([0.5, d/2 - 0.5])
    while coordenada[0] <= coordenada[1]:

        if thin == False:
            if d%2 != 0:
                posibilidades = [[coordenada[0], coordenada[1] - 1],[coordenada[0] + 1, coordenada[1]]]
                dists = [distancia_centro(posibilidades[0]), distancia_centro(posibilidades[1])]
                diferencia1 = np.abs(dists[0] - radio)
                diferencia2 = np.abs(dists[1] - radio)
            else:
                posibilidades = [[coordenada[0], coordenada[1] - 1],[coordenada[0] + 1, coordenada[1]]]
                dists = [distancia_centro(posibilidades[0]), distancia_centro(posibilidades[1])]
                diferencia1 = np.abs(dists[0] - radio)
                diferencia2 = np.abs(dists[1] - radio)

        elif thin == True:
            if d%2 != 0:
                posibilidades = [[coordenada[0] + 1, coordenada[1] - 1], [coordenada[0] + 1, coordenada[1]]]
                dists = [distancia_centro(posibilidades[0]), distancia_centro(posibilidades[1])]
                diferencia1 = np.abs(dists[0] - radio)
                diferencia2 = np.abs(dists[1] - radio)
            else:
                posibilidades = [[coordenada[0] + 1, coordenada[1] - 1],[coordenada[0] + 1, coordenada[1]]]
                dists = [distancia_centro(posibilidades[0]), distancia_centro(posibilidades[1])]
                diferencia1 = np.abs(dists[0] - radio)
                diferencia2 = np.abs(dists[1] - radio)
        
        temp = []
        if dists[0] > radio:
            if thin == False:
                temp.append(posibilidades[1][0])
                temp.append(posibilidades[1][1])
                coordenada[0] = posibilidades[1][0]
                coordenada[1] = posibilidades[1][1]
            
        elif dists[1] > radio:
            temp.append(posibilidades[0][0])
            temp.append(posibilidades[0][1])
            coordenada[0] = posibilidades[0][0]
            coordenada[1] = posibilidades[0][1]
            
        else:
            if dists[1] > dists[0]:
                temp.append(posibilidades[1][0])
                temp.append(posibilidades[1][1])
                coordenada[0] = posibilidades[1][0]
                coordenada[1] = posibilidades[1][1]
                
            else:
                temp.append(posibilidades[0][0])
                temp.append(posibilidades[0][1])
                coordenada[0] = posibilidades[0][0]
                coordenada[1] = posibilidades[0][1]
                
        bloques.append(temp)
    return bloques
         
def rotation(bloq, degrees):
    """
    This function will get a list of positions and then return the positions of all the rotated vectors.
    To perform the rotation we used the standard matrix form vector rotations in 2D.

    Input:
    ------
    bloq: list of lists representing points in the form [[x1, y1], [x2, y2], ..., [xn, yn]].
    degrees: float which represents the angle by which we want to rotate the points.

    Return:
    -------
    new_bloques: list of lists in the form [[x1, y1], [x2, y2], ..., [xn, yn]] containing the rotated
    points.

    """

    new_bloques = []
    for element in bloq:
        temp = []
        temp.append(element[0]*np.cos(degrees) - element[1]*np.sin(degrees))
        temp.append(element[0]*np.sin(degrees) + element[1]*np.cos(degrees))
        new_bloques.append(temp)

    return new_bloques

def is_circle(lista):
    """
    This function recibes a list of 2d points and returns a list of the distances of the points to the
    origin of the reference frame.

    Input:
    ------
    lista: list of lists of the form [[x1, y1], [x2, y2], ..., [xn, yn]] containing the points.
    
    Return:
    -------
    distancias: list of floats representing the distance of each of the blocks.

    """
    distancias = []
    for element in lista:
        dist = distancia_centro(element)
        distancias.append(dist)
    return distancias

def reflection(lista):
    """
    This function recives a list in the form [[x1, y1], [x2, y2], ..., [xn, yn]] and
    then returns a list in the same form but whith each point now reflected around the
    y axis.

    Input:
    ------
    lista: list of lists in the form [[x1, y1], [x2, y2], ..., [xn, yn]] containing the
    points to be rotated.

    Return:
    -------
    lista_nueva: list of lists in the form [[x1, y1], [x2, y2], ..., [xn, yn]] containing
    the reflected points.

    """
    lista_nueva = []
    
    for element in lista:
        temp = []
        temp.append(-element[0])
        temp.append(element[1])
        lista_nueva.append(temp)
    return lista_nueva

def main():
    """
    This is the main function of the programme.

    Input:
    ------
    None

    Return:
    -------
    None
    """
    # Preparacion
    radio, coordenada = radius_treatment()
    # Fin Preparacion

    # Creacion del 1/8 de circulos
    bloques = midpoint_algorithm(coordenada, radio)
    # Fin creacion del 1/8 de circulo

    bloques_totales = []
    bloques_totales_temp = []
 
    for element in bloques:
        bloques_totales.append(element)
        bloques_totales_temp.append(element)
    
    for elem in reflection(bloques):
        bloques_totales.append(elem)
        bloques_totales_temp.append(elem)


    for i in range(1,4): 
        temp = rotation(bloques_totales_temp,(np.pi/2) * i)
        for element in temp:
            bloques_totales.append(element)

    cantidad_neta = len(bloques_totales)
    if d % 2 == 0:
        cantidad = (cantidad_neta - 12) #la lista tiene 8 bloques repetidos y por lo rtanto debemos restarle 8 a su len(bloques_totales)
    else:
        cantidad_si_thick = cantidad_neta - 8
        cantidad = cantidad_neta - 12

    plot_list(bloques_totales)
    print("La cantidad de bloques son:", cantidad, "y la cantidad (No corregida) neta fue:" , cantidad_neta)
    print("The average radius is:", np.average(is_circle(bloques_totales)))
    print("The standard deviation of the radius is:",np.std(is_circle(bloques_totales)))
    print("The expected radius is:" , d/2)
    print("The measured number pi is:", (cantidad) / d)

    return None

main()
