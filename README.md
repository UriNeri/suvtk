# gb_submitter

[![PyPI](https://img.shields.io/pypi/v/gb_submitter.svg)](https://pypi.org/project/gb_submitter/)
[![Changelog](https://img.shields.io/github/v/release/LanderDC/gb_submitter?include_prereleases&label=changelog)](https://github.com/LanderDC/gb_submitter/releases)
[![Tests](https://github.com/LanderDC/gb_submitter/actions/workflows/test.yml/badge.svg)](https://github.com/LanderDC/gb_submitter/actions/workflows/test.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/LanderDC/gb_submitter/blob/master/LICENSE)

Tool to submit viral sequences to Genbank.

## Installation

Install this tool using `pip`:
```bash
pip install gb_submitter
```
## Usage

For help, run:
```bash
gb_submitter --help
```
You can also use:
```bash
python -m gb_submitter --help
```
## Development

To contribute to this tool, first checkout the code. Then create a new virtual environment:
```bash
cd gb_submitter
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
