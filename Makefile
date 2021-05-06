# --------------------------------------------------------- #
# *** WARNING ***: to use this Makefile make sure to change #
# the LIBS value to the path of your BOOST library!         #
# --------------------------------------------------------- #

CC = g++
CFLAGS = -Wall -O3 -Ofast --std=c++11
LIBS = /usr/lib/boost_1_76_0

default: main.exe

main.exe: main.o functions.o
	$(CC) $^ -o $@ -I $(LIBS)

%.o: %.cpp %.h
	$(CC) -c $< $(CFLAGS) -I $(LIBS)

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) -I $(LIBS)

clean :
	rm *.o *.exe