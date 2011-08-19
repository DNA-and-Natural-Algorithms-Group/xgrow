X11_FLAGS=-I/usr/X11R6/lib/include -L/usr/X11R6/lib -lX11

#if pkg-config is not in your default path, add the full path here
PKG_CONFIG=pkg-config

GLIB_CFLAGS=`${PKG_CONFIG} --cflags glib-2.0`

# For Mac OS / fink (may need to change /sw to /opt)
#GLIB_CFLAGS=-I/sw/lib -I/sw/include/glib-2.0 -I/sw/lib/glib-2.0/include


# For linux
GLIB_LIBS=`${PKG_CONFIG} --libs glib-2.0`

# For Mac OS / fink (may need to change /sw to /opt)
#GLIB_LIBS=-L/sw/lib -lglib-2.0 -lintl

xgrow: xgrow.c grow.c grow.h Makefile
	gcc -Wall -g -O3 -o  xgrow xgrow.c grow.c  ${X11_FLAGS} -lm 

xgrow-debug: xgrow.c grow.c grow.h Makefile
	gcc -Wall -g -o  xgrow xgrow.c grow.c  ${X11_FLAGS} -lm 

xgrow-small: xgrow.c grow.c grow.h Makefile
	gcc -Wall -g -O3 -o  xgrow-small xgrow.c grow.c -DSMALL ${X11_FLAGS} -lm 

xgrow-test: xgrow.c grow.c grow.h xgrow-tests.c xgrow-tests.h Makefile
	gcc -Wall  -O3 -g  -o  xgrow xgrow.c grow.c xgrow-tests.c -DTESTING_OK ${X11_FLAGS}  ${GLIB_CFLAGS} ${GLIB_LIBS} -lm 

clean: 
	rm -f xgrow xgrow-small 
