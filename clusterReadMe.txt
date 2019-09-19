!!! IN YOUR root ~ Directory !!!

> wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
> chmod +x Anaconda3-2019.07-Linux-x86_64.sh
> ./Anaconda3-2019.07-Linux-x87_64.sh

This is python 3.7
Hold enter through Terms of Service and type 'yes' to agree to them. Unless you don't agree, then game over
============================
Check if anaconda is on PATH
============================
> echo $PATH

If it is, you should see something like this somewhere in the output

/nas/longleaf/home/dsprague/anaconda3/bin

I'm 99% sure that it is added to your path by installation

IF NOT, do this, still in ~ (root) directory:

nano .bashrc

There should be a bunch of garbage in this file (not empty). At the bottom the file, add:

export PATH="/nas/longleaf/home/dsprague/anaconda3/bin:$PATH"

and save the file (control + x, yes yes)

============================
Check python version
============================

> which python

Should read ~/anaconda3/bin/python

> python -V

Should read python 3.7.X :: Anaconda custom (64-bit)

While python 3.7 should now be default, if you are like me you can ensure you are running python 3.7 by typing python3.7 instead of just python
Like so:

> python3.7 setup.py build_ext --inplace

============================
Standard mSEEKR install
============================
1. Dependencies
  a. cython
    Try "which cython", if none

    pip install cython

    OR

    conda install -c anaconda cython
  b. SEEKR
    pip install seekr

2. Type following commands in desired directory

  git clone https://github.com/spragud2/mSEEKR.git
  cd mSEEKR/
  python3.7 setup.py build_ext --inplace

3. Ignore warnings (unimportant GCC compiler warnings)
