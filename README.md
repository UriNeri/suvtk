# suvtk (Submission of Uncultivated Viral genomes toolkit)

[![PyPI](https://img.shields.io/pypi/v/suvtk.svg)](https://pypi.org/project/suvtk/)
[![Changelog](https://img.shields.io/github/v/release/LanderDC/suvtk?include_prereleases&label=changelog)](https://github.com/LanderDC/suvtk/releases)
[![Tests](https://github.com/LanderDC/suvtk/actions/workflows/test.yml/badge.svg)](https://github.com/LanderDC/suvtk/actions/workflows/test.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/LanderDC/suvtk/blob/master/LICENSE)

Tool to submit viral sequences to Genbank.

## Documentation

Documentation for the tool (including installation instruction) is available <a href="https://landerdc.github.io/suvtk/" target="_blank">here</a>.

## Usage

For help, run:
```bash
suvtk --help
```
You can also use:
```bash
python -m suvtk --help
```
## Development

To contribute to this tool, first checkout the code. Then create a new virtual environment:
```bash
cd suvtk
python -m venv venv
source venv/bin/activate
```
Now install the dependencies and test dependencies:
```bash
pip install -e '.[test]'
```
To run the tests:
```bash
python -m pytest
```
