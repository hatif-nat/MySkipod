# -*- coding: UTF-8 -*-
from sys import argv
import datetime
import os

key = argv[1][1:]

pows = range(25,31)
cores = [1, 2, 4, 8, 16, 32, 64, 128, 192, 256, 320, 384, 448, 512]

if(not os.path.exists('.tmp')):
    os.mkdir('.tmp')

tmp = '.tmp/' + key + '/'
if(not os.path.exists(tmp)):
    os.mkdir(tmp)

for pow in pows:
  print(pow)
  for core in cores:
     for times in range(0,3):
           if key == '-mpi':
                output_filename = tmp + 'mpi-'+str(core)+'-'+str(pow)+'.csv'
                if len(argv) == 3:
                    if argv[2] == '-pc':
                        system_call = 'mpirun -np' + str(core) + './test_mpi' +
                               + str(pow) + ' ' + str(output_filename)
                    else: 
                        print('Error key')
	            else:
                    system_call = 'mpisubmit.bg -n ' + str(core) + ' -w 5:00 ./test_mpi -- ' +
                               + str(pow) + ' ' + str(output_filename)
           elif key == '-omp':
                output_filename = tmp + 'omp-'+str(core)+'-'+str(pow)+'.csv' 
                system_call = './test_omp '+ str(core) + ' ' + str(pow)+' '+ str(output_filename)
           else:
               print('Error key')
	    k = os.system(system_call)

