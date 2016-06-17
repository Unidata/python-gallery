# Notebook Gallery
This is a collection of Jupyter notebooks demonstrating the use of Python
for meteorology and atmospheric science. The goal is to highlight and demonstrate
many of the useful libraries that exist in this space.

Some of these notebooks have been published on the Unidata
[Developer Blog](https://www.unidata.ucar.edu/blogs/developer/). In
addition to viewing those notebooks on the blog, you can also view all notebooks
on [nbviewer](http://nbviewer.jupyter.org/github/unidata/notebook-gallery/tree/master/).

## Installation Instructions

The current list of dependencies for this collection of notebooks is:

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
git clone https://github.com/Unidata/notebook-gallery

cd notebook-gallery

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

## Running the jupyter notebook server
From the `notebook-gallery` directory and with the `gallery` environment active,
run:

```sh
jupyter notebook
```

This starts the webserver for running jupyter notebooks. You should see some output like the following:

```
[I 11:59:31.409 NotebookApp] Serving notebooks from local directory: /Users/rmay/repos/notebook-gallery
[I 11:59:31.409 NotebookApp] 0 active kernels
[I 11:59:31.409 NotebookApp] The Jupyter Notebook is running at: http://localhost:8888/
[I 11:59:31.409 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
```

This indicates that the server is ready to accept connections on your local machine at port 8888.
Here's a link for convenience: [http://localhost:8888/](http://localhost:8888/).
If port 8888 is not available (say if you're running more than one server, you
may get a few more messages and the server may end up accepting connections on
a different port. The messages in the terminal should tell you which port to try.

Once you open up the main notebook page, if you click on the notebooks/ directory, you should
see the full collection of notebooks. Clicking on any of the notebooks will open them for running
in an interactive python session. Clicking a cell and typing shift-enter will run the cell.
The Jupyter Notebook [documentation](https://jupyter-notebook-beginner-guide.readthedocs.org)
has more information on working with notebooks.
