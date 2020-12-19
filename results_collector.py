# -*- coding: UTF-8 -*-
import datetime
import os

def get_general_table(path):
    general_table = []
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

    results = open('results/'+ key + '_' + str(datetime.datetime.now().isoformat()) + '.csv', 'w')


    sorted_table = []

    for node in general_table:
        sorted_table.append([int(node[0]), int(node[1]), float(node[2])])

    sorted_table = sorted(sorted_table)

    for node in sorted_table:
        results.write(','.join([str(key) for key in node]) + '\n')

    results.close()

if os.path.exists('.tmp/omp'):
    general_table = get_general_table('.tmp/omp/')
    create_general_table_file(general_table, 'omp')

if os.path.exists('.tmp/mpi'):
    general_table = get_general_table('.tmp/mpi/')
    create_general_table_file(general_table, 'mpi')