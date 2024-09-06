from setuptools import setup
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.md').read_text()

setup(
    name = 'sumofradicals',
    version = '0.9',
    description = 'A Python module to handle linear combinations of radicals exactly.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    
    author = 'Sebastian Gössl',
    author_email = 'goessl@student.tugraz.at',
    license = 'MIT',
    
    py_modules = ['sumofradicals'],
    url = 'https://github.com/goessl/sumofradicals',
    python_requires = '>=3.12',
    install_requires = ['sympy'],
    
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics'
    ]
)
