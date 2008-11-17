X11_FLAGS=-I/usr/X11R6/lib/include -L/usr/X11R6/lib -lX11

# If the below options don't work, replace them with the output of 

GLIB_CFLAGS=`pkg-config --cflags glib-2.0`

# For Mac OS / fink (may need to change /sw to /opt)
#GLIB_CFLAGS=-I/sw/lib -I/sw/include/glib-2.0 -I/sw/lib/glib-2.0/include


# For linux
GLIB_LIBS=`pkg-config --libs glib-2.0`

# For Mac OS / fink (may need to change /sw to /opt)
#GLIB_LIBS=-L/sw/lib -lglib-2.0 -lintl

xgrow: xgrow.c grow.c grow.h Makefile
	gcc -Wall   -g  -o  xgrow xgrow.c grow.c  ${X11_FLAGS} -lm 


xgrow-test: xgrow.c grow.c grow.h xgrow-tests.c xgrow-tests.h Makefile
	gcc -Wall  -O3 -g  -o  xgrow xgrow.c grow.c xgrow-tests.c -DTESTING_OK ${X11_FLAGS}  ${GLIB_CFLAGS} ${GLIB_LIBS} -lm 

