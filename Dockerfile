FROM informaticsmatters/rdkit-python-debian:Release_2018_09_1
ADD requirements.txt requirements.txt
USER root
RUN pip install -r requirements.txt
RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y git procps
RUN git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb
RUN pip install /usr/local/mmpdb 
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
