# ff1d-ftsw-pub-ice

Jupyter notebooks (Interactive computing environements, ICE) accompanying `SWIIFT
v0.10: A numerical model of wave-induced sea ice breakup based on an energy
criterion`.

## Table of contents

- [Installation](#installation)
- [Use](#use)


## Installation

Your Python version needs to be at least 3.11.

The `requirements.txt` file can be used to setup a virtual environement:
```console
python -m venv .venv
source .venv/bin/activate
python -m pip -r requirements.txt
```

The Jupyter server can they be launched with ```jupyter lab```.

> [!NOTE]
> We use SWIIFT v0.7 in this work. The namespace used by SWIIFT < v0.10 was
> flexfrac1d.


## Use

We provide four notebooks:
* `FractureThreshold.ipynb` can be used to run the numerical experiments;
* `Results.ipynb` can be used to reproduce the analysis of Section 4.1;
* `Results_ensemble.ipynb` can be used to reproduce the analysis of Section 4.2;
* `Schema.ipynb` can be used to reproduce the figures illustrating Section 2.
