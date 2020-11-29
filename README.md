### Task 59  
#### Разработка параллельной версии программы для вычисления определенного интеграла с использованием метода Симпсона

#### Compile
`cd build`
`make`

#### Запуск заданий
##### MPI 
`python job_loader.py -mpi`

##### OMP
`python job_loader.py -omp` 

#### Сбор результатов
`python results_collector.py` 

#### Построение графиков
`python analyzer.py`
