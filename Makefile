
xgrow: xgrow.c grow.c grow.h
	gcc -O -Wall -g -o xgrow xgrow.c grow.c -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lm

xgrow-test: xgrow.c grow.c grow.h
	gcc -D TESTING -O -Wall -g -o xgrow-test xgrow.c grow.c -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lm
