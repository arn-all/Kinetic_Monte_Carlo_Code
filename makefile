run: clean
	bin/DislocationKMC --config test.config

clean: build
	rm -f *.vtk

build: 
	make --directory=bin