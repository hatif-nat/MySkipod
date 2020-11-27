# from mpl_toolkits import mplot3d
from sys import argv
import datetime
import os
# import matplotlib.pyplot as plt

script, key = argv
# pows = range(10,26)
pows = [5]
# cores = [1]+list(range(2, 9, 2))
cores = [8]
if(not os.path.exists('tmp')):
    os.mkdir('tmp')


for pow in pows:
    for core in cores:
        for times in range(0,3):
            if key == '-mpi':
                outputFileName = 'tmp/mpi_{core}_{pow}.scv'.format( 
                    core = core, pow = pow
                )

                systemCall = 'mpirun -np {core} test_mpi {pow} {outputFileName}'.format(
                    core = core, pow = pow, outputFileName = outputFileName
                )
            else:
                outputFileName = 'tmp/omp_{core}_{pow}.scv'.format( 
                    core = core, pow = pow
                )

                systemCall = './test_omp {core} {pow} {outputFileName}'.format(
                    core = core, pow = pow, outputFileName = outputFileName
                )

            k = os.system(systemCall)

general_table = []

for outputFile in sorted(os.listdir('tmp')):
    file = open('tmp/' + outputFile, 'r')
    lines = file.readlines()
    cores = lines[0].split(',')
    times = [float(line.split(',')[2]) for line in lines]
    md_time = sum(times)/3
    node = lines[0].split(',')[:-1]
    node.append(str(md_time))
    general_table.append(node) 
    file.close()

os.system('rm -rf tmp')  

if(not os.path.exists('gt_save')):
    os.mkdir('gt_save')
g_t_safe = open('gt_save/'+ key[1:] + '_g_t_' + str(datetime.datetime.now()) + '.scv', 'w')

for node in general_table:
    g_t_safe.write(','.join(node) + '\n')

g_t_safe.close()

# num_bars = len(general_table)

# fig = plt.figure();

# ax = plt.axes(projection= "3d")

# ax.set_xlabel('Cores')
# ax.set_ylabel('Pow of 2')
# ax.set_zlabel('Runtime (s)')

# x_pos = [int(node[0]) for node in general_table]
# y_pos = [int(node[1]) for node in general_table]
# z_pos = [0] * num_bars

# z_size = [float(node[2]) for node in general_table]
# x_size = [0.25 for i in z_size]
# y_size = [0.25 for i in z_size]

# ax.bar3d(x_pos, y_pos, z_pos, x_size, y_size, z_size)

# plt.show()

