FROM python:3.11.4-slim-bullseye

USER root
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y \
        git \
        procps && \
    pip install rdkit==2023.3.2 && \
    git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb && \
    pip install /usr/local/mmpdb

COPY requirements.txt ./
RUN pip install -r requirements.txt

ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
