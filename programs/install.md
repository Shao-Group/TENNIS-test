# Installation
All programs should be installed/softlinked in this directory (`programs`).
## TENNIS
TENNIS should be installed and softlinked in this directory (i.e. `programs/tennis.py`).

For example, 
```sh
cd programs

# install dependency PySAT
pip install python-sat[aiger,approxmc,cryptosat,pblib]

# install TENNIS
git clone https://github.com/Shao-Group/TENNIS

# softlink
ln -s TENNIS/src/tennis.py ./tennis.py
```
## chr mapping

The easiest way to install `cvbio` is via conda. Read more about [cvbio](https://github.com/clintval/cvbio#cvbio).

```sh
conda install cvbio
```
Then we can create a softlink for it.
```sh
ln -s $(which cvbio) ./cvbio
```