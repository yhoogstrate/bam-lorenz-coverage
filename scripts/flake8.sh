#!/bin/bash

source .venv/bin/activate && flake8 --ignore=E501 blc/ tests/
