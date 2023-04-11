CC=g++

all: main.cpp
	$(CC) main.cpp -o Stocks

clean:
	rm -f Stocks