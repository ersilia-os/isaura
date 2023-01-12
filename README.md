# The Isaura data store

A library to cache precalculated properties of biomedical entities on remote infrastructure (DynamoDB) and locally (sqlite3).

This repository provides an interface to the precalculated data available from the Ersilia Model Hub. At the moment, Isaura is focused on chemical descriptors.

## Quick start guide

### 1. Clone the repository

```bash
git clone https://github.com/ersilia-os/isaura.git
cd isaura
```
### 2. Create a conda environment and activate it

```bash
conda env create -f env.yaml
conda activate isaura
```

### 3. Clone the repository and install it with pip

```bash
pip install -e .
```

### 4. Use Isaura to store pre-calculations locally
Create an `IsauraLocalClient` to store pre calculations fetched from the Ersilia Model Hub locally

```python
from isaura.service.client import IsauraLocalClient

# This will initialize a local sqlite3 database at ~/.local/eos/isaura_local.db
local_client = IsauraLocalClient()

# Insert precalcs in bulk
local_client.insert([precalc for precalc in precalcs])
```

_If you are looking for the old isaura version to locally store precalculations in .h5 format, please use the release v0.1_

### 5. Fetch pre-calculations from online server
Once Isaura is installed, you can start using Isaura Clients to fetch pre calculations from the Ersilia remote server and store them in your local cache. First create an `IsauraRemoteClient` to fetch pre calculations from ersilia

```python
from isaura.service.client import IsauraRemoteClient

# Initialize the user client with the API url
# Find the url for Ersilia Precalc API [here]
remote_client = IsauraRemoteClient(url = [Ersilia Precalc API URL])

# Client returns a `ResponseBodySchema` object
resp = remote_client.get_all_precalcs()
precalcs = resp.items
```

Please look at sections below for more detailed examples and documentation of the programming API.

## Isaura clients

### Isaura remote client

Remote client interacts with remote cache to perform read functions. Find the programming API below.

```python
from isaura.service.client import IsauraRemoteClient

# Initialize the user client with the API url
# Find the url for Ersilia Precalc API [here]
remote_client = IsauraRemoteClient(url = [Ersilia Precalc API URL])

# Client returns a `ResponseBodySchema` object
resp = remote_client.get_all_precalcs(last_eval_key=None)
precalcs = resp.items

# If last_eval_key is not `None` then more data is available
# pass last_eval_key to `get_all_precalcs` function to get rest of the data
last_eval_key = resp.last_eval_key

# Get precalc by id
resp = remote_client.get_precalc_by_id(precalc_id="")
precalc = resp.items[0]

# Get precalcs by model id
resp = remote_client.get_precalcs_by_model_id(model_id="")
precalcs = resp.items

# Get precalc by input key
# A model id list is required
resp = remote_client.get_precalcs_by_input_key(model_id_list=[], input_key = "")
precalcs = resp.items
```

### Isaura local client

Local client interacts with local cache to perform read, write and delete functions.

```python
from isaura.routes.schemas.common import Precalc
from isaura.service.client import IsauraLocalClient

# Local database can be created at a custom path
# Defaults to ~/.local/eos/isaura_local.db
local_client = IsauraLocalClient(db_path=[path to db])

# Reset the local database
local_client.reset()

# Use the Precalc class to create precalc objects
precalc = Precalc(model_id = "model id", input_key = "input key", value = {"out" : "model output value"})

# Insert precalcs in bulk
local_client.insert([precalc])

# Delete precalc in bulk with precalc ids
local_client.delete([precalc.precalc_id])

precalcs = local_client.get_all_precalcs(page = 0, limit = 100)

# Get precalc by id
precalc = local_client.get_precalc_by_id(precalc_id="")[0]

# Get precalcs by model id
precalcs = local_client.get_precalcs_by_model_id(model_id="")

# Get precalc by input key
precalcs = local_client.get_precalcs_by_input_key(input_key = "")
```

### Isaura Admin client

Admin client interacts with remote cache to perform insert and delete functions. An AWS account with permissions to isaura dynamo table is required to use admin client.

> You can skip this section if you only want to fetch precalculations hosted by Ersilia.

```python

from isaura.routes.schemas.common import Precalc
from isaura.service.client import IsauraAdminClient

admin_client = IsauraAdminClient()

# Use the Precalc class to create precalc objects
precalc = Precalc(model_id = "model id", input_key = "input key", value = {"out" : "model output value"})

# Insert precalcs in bulk
admin_client.insert([precalc])

# Delete precalc in bulk with precalc ids
admin_client.delete([precalc.precalc_id])
```

## Provision AWS infrastructure

This section explains how to host your own aws infrastructure for remote cahce (Do you want to host your own DynamoDB and pay for it?).

**If you just want to use the pre calculations hosted by Ersilia then you can skip this section.**

### Requirements

- Python >=v3.7
- Nodejs >=v18.12.1

### Create a conda environment and activate it

```bash
conda env create -f env.yaml
conda activate isaura
```

### Create a python virtual env

`.venv_prod` is used for creating lambda layers used by lamda function on aws.
Do not use this virtual envirinment for anything else.

```bash
# ! This is important. Do not use conda or any other venv alternatives
python -m venv .venv_prod
```

### Install CDK CLI dependencies

```bash
npm install
```

### Install Python dependencies

```bash

# build isaura package for prod env
python -m build

# activate prod venv and install packaage
source .venv_prod/bin/activate
pip install dist/[lastest built wheel]
deactivate
# This is required so that our lambda functions can use isaura package
```

### Bootstrap CDK

This only needs to be done once for an aws account.

```bash
npx cdk bootstrap
```

### Deploy Isaura infrastructure

```bash
npx cdk deploy
```

## License

This repository is open-sourced under [the GPL-3 License](https://github.com/ersilia-os/isaura/blob/master/LICENSE).

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission or [volunteer](https://www.ersilia.io/volunteer) with us!
