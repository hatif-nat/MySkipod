# -*- coding: UTF-8 -*-
from sys import argv
import datetime
import os

key = argv[1][1:]

sizes = [128, 256, 512, 1024] #, 1536, 2048, 2560, 3072, 3584]
cores = [2, 4, 8, 12] #, 16, 32, 64, 128, 160]

if(not os.path.exists('.tmp')):
    os.mkdir('.tmp')

tmp = '.tmp/' + key + '/'
if(not os.path.exists(tmp)):
    os.mkdir(tmp)
sys_call = 'echo "Hi"'
for size in sizes:
    for core in cores:
        for times in range(0,3):
            if (key == 'mpi'):
                output_filename = tmp + 'mpi-' + str(core) +'-' + str(size) + '.csv'
                if (len(argv) == 3):
                    if argv[2] == '-pc':
                        sys_call = 'mpirun -np ' + str(core) + ' ./test_mpi ' + str(size) + ' ' + str(output_filename)
                    else: 
                        print('Error key')
                else: 
                    sys_call = 'mpisubmit.pl -p ' + str(core) + ' -w 00:05 ./test_mpi -- ' + str(size) + ' ' + str(output_filename)
            elif (key == 'omp'):
                output_filename = tmp + 'omp-' + str(core) +'-' + str(size) + '.csv'
                sys_call = './test_omp ' + str(core) + ' ' + str(size) + ' ' + str(output_filename)
            else:
                print('Error key')
            print(sys_call)
            k = os.system(sys_call)

