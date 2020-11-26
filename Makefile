FILES=test.cpp
APP_NAME=t
FLAGS=-fopenmp -o

all:
	g++ $(FILES) $(FLAGS)   $(APP_NAME)

clean:
	rm -rf $(APP_NAME)