FILES=test.cpp
APP_NAME=t
OPTIMIZE_LEVEL=3
FLAGS=-Wall -Wunreachable-code -pedantic

all:
	g++ $(FLAGS) -O$(OPTIMIZE_LEVEL) $(FILES) -o $(APP_NAME)

clean:
	rm -rf $(APP_NAME)