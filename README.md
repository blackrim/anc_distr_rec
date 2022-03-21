# continuous character distribution reconstructions

# requirements
This is run with python 3 (because of some shifts in the libraries). This requires `numpy`, `scipy`, and `matplotlib`. You can add viz with `ete3`.

I installed these things like this

- `pip3 install numpy scipy matplotlib`
- `python3 -m pip install PyQt5`
- `pip3 install ete` (may have to do `pip3 install --user ete3`)

# running an example
You can run the example with `python3 src/main.py examples/test.cont.tre`. Right now, this will simulate the data with a 20% of a bimodal distribution and 80% of Rayleigh. To run your own data, it needs to be in phylip format (this is changing to csv in just a sec) and run `python src/main.py examples/test.cont.tre examples/test.cont`.

# contributors
 - Stephen A. Smith
 - Brian O'Meara

(Some of this has been incorporated into the biotaphy project)
