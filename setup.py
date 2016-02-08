from setuptools import setup, Extension
import numpy

config = {
    'description': 'xgrowutils',
    'author': 'Constantine Evans',
    'author_email': 'cgevans@evans.foundation',
    'version': '0.2.0',
    'install_requires': ['numpy', 'pyyaml'],
    'extras_require': { 'parallel': 'ipyparallel' },
    'packages': ['xgrowutils', 'xgrowutils.scripts'],
    'entry_points': { 'console_scripts':
        ['xgrow-wrap = xgrowutils.scripts.xgrow_wrap:main']
        },
    'name': 'xgrowutils'
}


setup(test_suite='nose.collector', **config)
