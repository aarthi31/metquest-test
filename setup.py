from setuptools import setup, find_packages

requirements = [
    'cobra>=0.5.4',
    'numpy',
    'scipy',
    'python-libsbml',
    'networkx>=2.1',
    'click'    
]

try:
    with open('README.md') as handle:
        description = handle.read()
except:
    description = ''

setup(
    name='metquest',
    version='0.1.0',
    packages=find_packages(),
    setup_requires=[],
    install_requires=requirements,
    tests_require=['pytest'],
    author='Aarthi Ravikrishnan', 'Meghana Nasre', 'Karthik Raman',
    author_email='aarthiravikrishnan@gmail.com',
    description='MetQuest: Enumerating all possible biosynthetic pathways in metabolic networks ',
    long_description=description,
    license="LGPL/GPL v2+",
    keywords='metabolism biology graph-theory pathways',
    url='https://github.com/RamanLab/MetQuest',
    download_url='https://testpypi.python.org/pypi/metquest',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2'
            ' or later (LGPLv2+)',
        'License :: OSI Approved :: GNU General Public License v2'
            ' or later (GPLv2+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms='GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7'
)

