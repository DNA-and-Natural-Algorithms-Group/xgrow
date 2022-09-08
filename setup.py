#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
from distutils.command.build import build
from setuptools.command.develop import develop

#import sys
#sys.argv.extend(['plat-name', 'x86_64'])

BUILD_STRING = (
    "{} -Wall -Wno-unused-result -g -O2 src/xgrow.c src/grow.c -o xgrow/_xgrow -lm {}"
)
BUILD_STRING_LIB = "{} -Wall -Wno-unused-result -g -O2 src/grow.c -fPIC -shared -o xgrow/libxgrow.so -lm"


def find_x11():
    import os

    if "X11_FLAGS" in os.environ:
        x11s = os.environ["X11_FLAGS"]
    else:
        x11f = []
        includes = ["/usr/include/X11", "/opt/X11/include"]
        for x in includes:
            # if x is None:
            #    raise Exception("Can't find an X11 include dir.")
            if os.path.exists(x):
                x11f.append("-I{}".format(x))
                break
        libs = ["/usr/lib/X11", "/opt/X11/lib", "/usr/lib64", "/usr/lib"]
        for x in libs:
            # if x is None:
            #    raise Exception("Can't find an X11 lib dir.")
            if os.path.exists(x):
                x11f.append("-L{}".format(x))
                break
        x11f.append("-lX11")
        x11s = " ".join(x11f)
    if "CC" in os.environ:
        cc = os.environ["CC"]
    else:
        cc = "cc"
    return (cc, x11s)


class build_xgrow(build):
    def run(self):
        import os

        os.system(BUILD_STRING.format(*find_x11()))
        os.system(BUILD_STRING_LIB.format(find_x11()[0]))

        build.run(self)


class develop_xgrow(develop):
    def run(self):
        import os

        os.system(BUILD_STRING.format(*find_x11()))
        os.system(BUILD_STRING_LIB.format(find_x11()[0]))

        develop.run(self)


setup(
    name="xgrow",
    version="20220725",
    packages=["xgrow"],
    install_requires=["pyyaml", "typing_extensions", "numpy", "pandas", "matplotlib"],
    include_package_data=True,
    package_data={"xgrow": ["_xgrow", "_xgrow.so", "py.typed"]},
    cmdclass={"build": build_xgrow, "develop": develop_xgrow},
    entry_points={"console_scripts": ["xgrow = xgrow._script_xgrow:main"]},
    author="Constantine Evans et al (this version)",
    author_email="const@costi.eu",
    description="Xgrow in pythonic form",
    zip_safe=False,
    has_ext_modules=lambda: True
)
