
# Compiler
CC = g++
# Compiler flags
CFLAGS = -lgsl -lgslcblas -std=c++17 -lboost_program_options

# Targets and object files for g-plot and test
G_PLOT_TARGET = g-plot
G_PLOT_OBJS = g-plot.o integral.o flag_parser.o error_handler.o variables.o main_functions.o

TEST_TARGET = test
TEST_OBJS = test.o integral.o flag_parser.o error_handler.o variables.o main_functions.o

all: $(G_PLOT_TARGET)

$(G_PLOT_TARGET): $(G_PLOT_OBJS)
	$(CC) $(CFLAGS) -o $(G_PLOT_TARGET) $(G_PLOT_OBJS)

$(TEST_TARGET): $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $(TEST_TARGET) $(TEST_OBJS)

test.o: test.cpp ./src/integral.h
	$(CC) $(CFLAGS) -c test.cpp

g-plot.o: g-plot.cpp ./src/integral.h
	$(CC) $(CFLAGS) -c g-plot.cpp

integral.o: ./src/integral.cpp ./src/integral.h
	$(CC) $(CFLAGS) -c ./src/integral.cpp

flag_parser.o: ./src/flag_parser.cpp ./src/flag_parser.h
	$(CC) $(CFLAGS) -c ./src/flag_parser.cpp

error_handler.o: ./src/error_handler.cpp ./src/error_handler.h
	$(CC) $(CFLAGS) -c ./src/error_handler.cpp

variables.o: ./src/variables.cpp ./src/variables.h
	$(CC) $(CFLAGS) -c ./src/variables.cpp

main_functions.o: ./src/main_functions.cpp ./src/main_functions.h
	$(CC) $(CFLAGS) -c ./src/main_functions.cpp

clean:
	rm -f *.o $(G_PLOT_TARGET) $(TEST_TARGET)
