all: movie

DIR := $(shell date +%T-%F)

movie: 
	convert -adjoin -loop 0 -delay 20 impstep*.png impstepMOVIE.gif
	mkdir $(DIR)
	mv impstep*.* $(DIR)/

clean:
	rm *.png
