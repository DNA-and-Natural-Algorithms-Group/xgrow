
xgrow: xgrow.c grow.c grow.h
	gcc -O -Wall -g -o xgrow xgrow.c grow.c -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lm
