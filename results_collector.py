# -*- coding: UTF-8 -*-
import datetime
import os

if os.path.exists('.tmp/omp')):
    general_table_omp = get_general_table('./tmp/omp'))
    create_general_table_file(general_table, 'omp')

if os.path.exists('.tmp/mpi')):
    general_table_omp = get_general_table('./tmp/mpi'))
    create_general_table_file(general_table, 'mpi')

def get_general_table(path):
    for outputFile in sorted(os.listdir(path)):
        file = open( path + outputFile, 'r')
        lines = file.readlines()
        times = [float(line.split(',')[2]) for line in lines]
        md_time = sum(times) / len(times)
        node = lines[0].split(',')[:-1]
        node.append(str(md_time))
        general_table.append(node) 
        file.close()
    return general_table

def create_general_table_file(general_table, key):
    if(not os.path.exists('results')):
        os.mkdir('results')
    
    results = open('results/'+ key + str(datetime.datetime.now()) + '.csv', 'w')

    general_table = sorted[ [(int)node[0], (int)node[1], node[2]] for node in general_table]

    for node in general_table:
        results.write(','.join(node) + '\n')

    results.close()
