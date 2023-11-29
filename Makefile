CC = mpicxx
SRCS = ./src/*.cpp
INC = -I ./src
OPTS = -std=c++17 -Wall -Werror -lpthread -O3
LINK = -lglfw3 -lGLEW -ldl -lpthread -lm -lGL -lGLU -lX11  
EXEC = bin/nbody
LIB = -L /usr/local/lib/libglfw3.a
DEF ?= NA
all: clean compile

compile:
	$(CC) -D${DEF} $(SRCS) $(OPTS) $(INC) $(LINK) -o $(EXEC) 

clean:
	rm -f $(EXEC)
