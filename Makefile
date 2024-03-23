TARGET = ./bin/pcg
CC = gcc 
CFLAGS = -O3 -I./include -I./libs/fft/include -I./libs/tri-dpf/include -I/opt/homebrew/opt/openssl/include
LDFLAGS = -march=native -lcrypto -lssl -lm -maes -ffast-math

# Define source files and filter out library files 
FFT_SRC = $(filter-out ./libs/fft/src/test.c, $(wildcard ./libs/fft/src/*.c))
DPF_SRC = $(filter-out ./libs/tri-dpf/src/test.c, $(wildcard ./libs/tri-dpf/src/*.c))
PCG_SRC = $(wildcard ./src/*.c)

all: $(TARGET)

OBJECTS = $(FFT_SRC:.c=.o) $(DPF_SRC:.c=.o) $(PCG_SRC:.c=.o)

$(TARGET): $(OBJECTS)
	@mkdir -p ./bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f libs/fft/src/*.o libs/tri-dpf/src/*.o *.o $(TARGET) $(OBJECTS)

.PHONY: all clean