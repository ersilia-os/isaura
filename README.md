# The Isaura data store

A library to cache precalculated properties of biomedical entities on remote infrastructure (DynamoDB) and locally (sqlite3).

This repository provides an interface to the precalculated data available from the Ersilia Model Hub. At the moment, Isaura is focused on chemical descriptors.

## Getting Started

### Requirement for local development

- Python >=v3.8
- Nodejs >=v18.12.1

### Create a conda environment

```bash
conda env create -f env.yaml
conda activate isaura
```

## Isaura clients

This library provide user friendly clients to interact with remote and local cache.

### Isaura Admin client

Admin client interacts with remote cache to perform insert and delete functions. An AWS account with permissions to isaura dynamo table is required to use admin client.

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

### Isaura user client

User client interacts with remote cache to perform read functions. No AWS account is required. read API is open for public use.

```python
from isaura.service.client import IsauraClient

# Initialize the user client with the API url
user_client = IsauraClient(url = "")

# Client returns a `ResponseBodySchema` object
resp = user_client.get_all_precalcs(last_eval_key=None)
precalcs = resp.items

# If last_eval_key is not `None` then more data is available
# pass last_eval_key to `get_all_precalcs` function to get rest of the data
last_eval_key = resp.last_eval_key

# Get precalc by id
resp = user_client.get_precalc_by_id(precalc_id="")
precalc = resp.items[0]

# Get precalcs by model id
resp = user_client.get_precalcs_by_model_id(model_id="")
precalcs = resp.items

# Get precalc by input key
# A model id list is required
resp = user_client.get_precalcs_by_input_key(model_id_list=[], input_key = "")
precalcs = resp.items
```

### Isaura local client

Local client interacts with local cache to perform read, write and delete functions.

```python
from isaura.routes.schemas.common import Precalc
from isaura.service.client import IsauraLocalClient

local_client = IsauraLocalClient()

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

## Provision AWS infrastructure

This is only required if you want to host you own remote cache.

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
