#!/usr/bin/env python

from setuptools import setup, find_packages
from distutils.command.build import build
from setuptools.command.develop import develop

BUILD_STRING = "{} -Wall -g -O2 src/xgrow.c src/grow.c -o xgrow/_xgrow -lm {}"

def find_x11():
    import os
    if 'X11_FLAGS' in os.environ:
        x11s = os.environ['X11_FLAGS']
    else:
        x11f = []
        includes = ['/usr/include/X11','/opt/X11/include',None]
        for x in includes:
            if x is None:
                raise Exception("Can't find an X11 include dir.")
            if os.path.exists(x):
                x11f.append("-I{}".format(x))
                break
        libs = ['/usr/lib/X11','/opt/X11/lib','/usr/lib64','/usr/lib',None]
        for x in libs:
            if x is None:
                raise Exception("Can't find an X11 lib dir.")
            if os.path.exists(x):
                x11f.append("-L{}".format(x))
                break
        x11f.append('-lX11')
        x11s = " ".join(x11f)
    if 'CC' in os.environ:
        cc = os.environ['CC']
    else:
        cc = 'cc'
    return (cc,x11s)
        
class build_xgrow(build):
    def run(self):
        import os
        os.system(BUILD_STRING.format(*find_x11()))
        
        build.run(self)

class develop_xgrow(develop):
    def run(self):
        import os
        os.system(BUILD_STRING.format(*find_x11()))
        
        develop.run(self)

setup(
    name = "xgrow",
    version = "20170507dev0",
    packages = ['xgrow'],

    install_requires = [],

    include_package_data=True,
    package_data= {'xgrow': ['_xgrow']},

    cmdclass={'build': build_xgrow, 'develop': develop_xgrow},
    
    entry_points={ 'console_scripts': [
        'xgrow = xgrow._script_xgrow:main']},

    author = "Constantine Evans et al (this version)",
    author_email = "cgevans@evans.foundation",
    description = "Xgrow in pythonic form"
)
