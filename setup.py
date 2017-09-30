from setuptools import setup, find_packages
from os.path import join, dirname

setup(
    name='ocelot',
    version='17.07rc',
    description='Accelerator, radiation and x-ray optics simulation framework',
    author='ocelot-collab',
    author_email='tomin.sergey@gmail.com',
    url='https://github.com/ocelot-collab/ocelot',
    packages=find_packages(),
    # long_description=open(join(dirname(__file__), 'README.txt')).read(),
)