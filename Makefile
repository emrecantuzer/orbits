orbits: orbits.o orbits.h
	gcc -o orbits orbits.o -lm -lX11 -lplplotd

orbits.o: orbits.c orbits.h
	gcc -c orbits.c
