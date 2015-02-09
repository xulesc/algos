CC=g++
CFLAGS=--shared -fPIC -O3 -ffast-math  -I/usr/include/python2.7/ -lpython2.7 -lm

munkres.so: munkres.c
	$(CC) -D __PYMOD__ -o $@ $< $(CFLAGS)

run: munkers.so psalign.py
	python psalign.py

clean:
	rm munkres.so
