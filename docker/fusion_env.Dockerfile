# Minimal variant for just postprocessing / merge
FROM python:3.10-slim

RUN pip install --no-cache-dir pandas numpy

WORKDIR /workspace
