# The Isaura data store

A lake of precalculated properties of biomedical entities

This repository provides an interface to the precalculated data available from the Ersilia Model Hub. At the moment, Isaura is focused on chemical descriptors.

## Getting Started

### Requirement for local development

- Python >=v3.8
- Nodejs >=v18.12.1

### Create two python virtual env

`.venv_dev` is used for development and local usage.

`.venv_prod` is used for creating lambda layers for production use.

```bash
# ! This is important. Do not use conda or any other venv alternatives
python -m venv .venv_dev
python -m venv .venv_prod
```

### Install poetry

Do this for both environments

```bash
pip install poetry
```

### Install CDK CLI dependencies

```bash
npm install
```

### Install Python dependencies

Install

```bash
# In dev environment
poetry install

# In production environment
poetry install --no-dev
```
