#!/bin/bash
set -e

uv run python cleanup.py
quarto render
git add -A
git commit -m "."
git push origin main
