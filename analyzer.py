# -*- coding: UTF-8 -*-
from mpl_toolkits import mplot3d
import os
import matplotlib.pyplot as plt


for file in sorted(os.listdir('results'))
    show_graph(filename)

def show_graph(filename):
    file = open(filename)

    lines = file.readlines()
    general_table = []

    for line in lines:
        general_table.append(line.split(','))

    num_bars = len(general_table)

    fig = plt.figure();

    ax = plt.axes(projection= "3d")
    ax.set_name(filename)
    ax.set_xlabel('Cores/Threads')
    ax.set_ylabel('Pow of 2')
    ax.set_zlabel('Runtime (s)')

    x_pos = [int(node[0]) for node in general_table]
    y_pos = [int(node[1]) for node in general_table]
    z_pos = [0] * num_bars

    z_size = [float(node[2]) for node in general_table]
    x_size = [10 for i in z_size]
    y_size = [0.7 for i in z_size]

    ax.bar3d(x_pos, y_pos, z_pos, x_size, y_size, z_size)

    plt.show()

