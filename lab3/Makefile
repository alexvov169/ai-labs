SOURCES := main.cpp hopfield_NN.cpp
BIN := main
# Objs are all the sources, with .cpp replaced by .o
OBJS := $(SOURCES:.cpp=.o)

all: $(BIN)

check: $(BIN)
	./$(BIN)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BIN): $(OBJS)
	g++ $(CFLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(BIN) $(OBJS)
