# queries
SPARQL queries for use with grlc

## Installation instructions
Pull requests on CLARIA/grlc are pending, until these are merged, a it is
necessary to use c-martinez/grlc fork. Run the following commands to create
a virtual environment and install grlc inside it:

```
virtualenv venv
source venv/bin/activate
pip install git+https://github.com/c-martinez/grlc.git@candyGene
```

You will need a `config.ini` file in your current working folder. Download it from [here](https://github.com/c-martinez/grlc/blob/candyGene/config.default.ini) and edit it
as required.

 <edit config.ini>

Grlc requires an external npm package (git2prov) to generate provenance information.
Install it like this:


```
npm install git2prov
```

Now you should be ready to run a grlc server locally. Start it with:

```
grlc-server
```

Browse to: `http://localhost:8088/api/candYgene/queries` and you should get a
local copy of the grlc candyGene API.

# Dockerized version

Instructions:
 - Clone this repo.

 - Run:
```
docker build -t candygene .
docker run -it --rm -p 8888:8888 -p 8088:8088 candygene
```

 - Browse to: `http://localhost:8088/api/local` and you should get a local copy of queries in this repo.

 - Browse to: `http://localhost:8888/` and you should see your notebook.
