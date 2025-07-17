from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name='ocelot-collab',
    version='25.07.0',
    description='Accelerator, radiation and x-ray optics simulation framework',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='ocelot-collab',
    author_email='sergey.tomin@desy.de',
    url='https://github.com/ocelot-collab/ocelot',
    packages=find_packages(include=["ocelot", "ocelot.*", "demos", "demos.*"]),
    package_dir={
        'ocelot': 'ocelot',
        'demos': 'demos'
    },
    install_requires=[
        'numpy', 'scipy', 'matplotlib', 'pandas', 'h5py'
    ],
    extras_require={'docs': ['Sphinx', 'alabaster', 'sphinxcontrib-jsmath']},
    package_data={
        'ocelot.optics': ['data/*.dat'],
        'ocelot': ['py.typed']
    },
    license="GNU General Public License v3.0",
    python_requires=">=3.9"
)
