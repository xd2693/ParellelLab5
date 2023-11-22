CC = g++ -std=c++11
SRCS = ./src/*.cpp
INC = -I ./src
OPTS = -std=c++17 -Wall -Werror -lpthread -O3
LINK = -lglfw3 -lGLEW -ldl -lpthread -lm -lGL -lGLU -lX11  
EXEC = bin/MPI
LIB = -L /usr/local/lib/libglfw3.a

all: clean compile

compile:
	$(CC) $(SRCS) $(OPTS) $(INC) $(LINK) -o $(EXEC)

clean:
	rm -f $(EXEC)
