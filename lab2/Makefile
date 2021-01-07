SOURCES := main.cpp aco.cpp ant.cpp
BIN := main
# Objs are all the sources, with .cpp replaced by .o
OBJS := $(SOURCES:.cpp=.o)

all: $(BIN)

check: $(BIN)
	printf "runall\n" | ./$(BIN)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BIN): $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(BIN) $(OBJS)
