# Unidata Python Gallery
This is a collection of examples published in the Unidata
[blogs](https://www.unidata.ucar.edu/blogs/developer/) or just contributed. In
addition to viewing them on the blog you can view them rendered at the [web gallery](http://unidata.github.io/python-gallery),
where you can also download the examples as Jupyter notebooks.

## Installation Instructions

The current list of dependencies for this collection is:

- [numpy](http://www.numpy.org/)
- [netcdf4-python](https://unidata.github.io/netcdf4-python/)
- [matplotlib](http://matplotlib.org/)
- [cartopy](http://scitools.org.uk/cartopy/)
- [siphon](http://siphon.readthedocs.org)
- [MetPy](http://metpy.readthedocs.org)

The easiest way to install these libraries is with [conda](http://conda.pydata.org/).

1. [Install Miniconda (Python 3.4) from Continuum Analytics](http://conda.pydata.org/miniconda.html).
  ([Determine if your OS 32 or 64 bit](http://www.akaipro.com/kb/article/1616#os_32_or_64_bit))
2. Once Miniconda is installed, from the command line (e.g., OS X terminal,
  cmd.exe), run these instructions to clone the repository and create the environment:

```sh
git clone https://github.com/Unidata/python-gallery

cd python-gallery

conda env create -f environment.yml
```

### From a Unix command line (e.g., OS X terminal)
If your default shell is NOT bash, first type `bash`.
To activate or switch to a conda environment, you can `source activate
<environment>`. For example,

```sh
source activate gallery
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
activate gallery
```

To switch and/or deactivate environments:

```sh
deactivate
activate <environment>
```
