#!/bin/bash

source .venv/bin/activate && flake8 --ignore=E501 --select=E,W,N blc/ tests/
