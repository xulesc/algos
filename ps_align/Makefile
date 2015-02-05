CC=g++
CFLAGS=--shared -fPIC -O3 -ffast-math  -I/usr/include/python2.7/ -lpython2.7 -lm

munkers.so: munkers.c
	$(CC) -D __PYMOD__ -o $@ $< $(CFLAGS)

run: munkers.so psalign.py
	python psalign.py

clean:
	rm munkers.so
