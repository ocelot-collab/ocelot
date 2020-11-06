from setuptools import setup, find_packages
from os.path import join, dirname

all_packages = []
for pkg in find_packages():
    if "demos" in pkg:
        pkg = "ocelot." + pkg
    all_packages.append(pkg)


setup(
    name='ocelot',
    version='20.11.2',
    description='Accelerator, radiation and x-ray optics simulation framework',
    author='ocelot-collab',
    author_email='tomin.sergey@gmail.com',
    url='https://github.com/ocelot-collab/ocelot',
    packages=all_packages,
    package_dir={'ocelot.demos': 'demos'},  ## install examples along with the rest of the source
    install_requires=[
        'numpy', 'scipy', 'matplotlib', 'pandas'
    ],
    extras_require={'docs': ['Sphinx', 'alabaster', 'sphinxcontrib-jsmath']},
    package_data={'ocelot.optics': ['data/*.dat']},
    license="GNU General Public License v3.0",
)
