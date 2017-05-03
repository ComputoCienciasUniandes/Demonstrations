from distutils.core import setup, Extension

module1 = Extension('astrohut/bruteforce',
                    sources = ['astrohut/bruteforce.c', 'astrohut/box.c', 'astrohut/init.c'], include_dirs=['astrohut'], extra_compile_args=["-fopenmp"],
                     extra_link_args=["-fopenmp"])

setup (name = 'AstroHut',
       version = '1.0',
       description = 'Bruteforce attemp',
       packages = ['astrohut'],
       ext_modules = [module1])
