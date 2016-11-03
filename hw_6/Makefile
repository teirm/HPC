all: particles particles_serial

particles: particles.c
	mpicc -std=c99  -o particles particles.c

particles_serial: particles_serial.c
	g++ -o particles_serial particles_serial.c

clean:
	rm -f particles
	rm -f particles_serial
