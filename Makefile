
xgrow: xgrow.c grow.c grow.h xgrow-tests.c xgrow-tests.h Makefile
	gcc -O -Wall -g -o xgrow xgrow.c grow.c xgrow-tests.c -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -lm -lglib-2.0 -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include

