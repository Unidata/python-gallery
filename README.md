# Blog Notebooks
This is a collection of notebooks and other resources published in the Unidata blogs.

## Installation Instructions

The current list of dependencies for this collection of notebooks is:

- [numpy](http://www.numpy.org/)
- [netcdf4-python](https://github.com/Unidata/netcdf4-python)
- [matplotlib](http://matplotlib.org/)
- [cartopy](http://scitools.org.uk/cartopy/)
- [siphon](https://github.com/Unidata/siphon)
- [MetPy](https://github.com/metpy/MetPy)

The easiest way to install these libraries is with [conda](http://conda.pydata.org/).

1. [Install Miniconda (Python 3.4) from Continuum Analytics](http://conda.pydata.org/miniconda.html).
  ([Determine if your OS 32 or 64 bit](http://www.akaipro.com/kb/article/1616#os_32_or_64_bit))
2. Once Miniconda is installed, from the command line (e.g., OS X terminal,
  cmd.exe), run these instructions to clone the repository and create the environment:

```sh
git clone https://github.com/Unidata/blog-notebooks

cd blog-notebooks

conda env create -f environment.yml
```

### From a Unix command line (e.g., OS X terminal)
If your default shell is NOT bash, first type =bash=.
To activate or switch to a conda environment, you can =source activate
<environment>=. For example,

```sh
source activate unidata-blog
```

To switch and/or deactivate environments:

```sh
source deactivate
source activate <environment>
```

### From a Windows command line (e.g., cmd.exe)

To activate or switch to a conda environment, you can `activate
<environment>`. For example,

```sh
activate unidata-blog
```

To switch and/or deactivate environments:

```sh
deactivate
activate <environment>
```
