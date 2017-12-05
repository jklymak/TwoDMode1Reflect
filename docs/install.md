---
layout: default
title: Installing
---
# Installing Software on a Mac

Tested on macOS 10, but should work on relatively recent OS X as well.

## Compiler code
   1. Install xcode using the App Store
   1. run `xcode-select --install` from command line
   1. Install [homebrew](http://brew.sh)
   3. install gcc (compilers) `brew install gcc --without-multilib`.
   4. install mpich (multi-core processing) `brew install mpich2`
   5. install hdf5 with mpi: `brew install hdf5 --with-mpi --with-fortran`
   6. intstal netcdf (file storage) `brew install -v netcdf --with-fortran`

Note that for the MITgcm as released `netcdf` does not require `hdf5` to be compiled with `mpi` (though the `MITgcm` itself requires `mpi`).  However, I am developing a parallel-aware netcdf interface for the `MITgcm`, so this step is useful.

## Python

   1. Install [anaconda](https://www.continuum.io/downloads)  I use python 2.7 still, but I think just about everything now works w/ python 3.5
   2. Try `ipython --pylab`.  From that command line try something like `fig,ax=plt.subplots();ax.plot(np.arange(10))`.
   3. Try `jupyter-notebook` and open the webpage they suggest.  See their [website](http://jupyter.org) for demos and I'll give some example python notebooks soon.
   4. install [xmitgcm](http://xmitgcm.readthedocs.io/en/latest/)
      1. `conda install xarray dask`
      2. `pip install xmitgcm`

## MITgcm

   1. Download the [source code](http://mitgcm.org/public/source_code.html).  Its a good idea to put this in a directory with the name the version you are running to avoid compatibility problems.  i.e. I keep a version at `~/MITgcmc66b/MITgcm` as that is the last version I downloaded.  I download the tar balls as they have been tested.  I always get errors from the CVS snapshots.
   2. Clone my example project

   ```
   git clone https://github.com/jklymak/MITgcmExampleSteadyGauss.git
   ```

  and then follow the instructions at the project [README](https://github.com/jklymak/MITgcmExampleSteadyGauss/blob/master/README.md)
