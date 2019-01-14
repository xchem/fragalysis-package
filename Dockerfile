FROM informaticsmatters/rdkit-python-debian:Release_2018_09_1
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
USER root
RUN apt-get update && apt-get install -y git procps
RUN git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb
RUN pip install /usr/local/mmpdb 
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis

# Conveneint command for built-in MolPort Neo4J processing...
#
#WORKDIR /usr/local/fragalysis/frag/network/scripts
#CMD ["./process_molport_compounds.py", \
#     "/exports/nextflow/fragbuilder/analysis/molport/2018-11", \
#     "iis_smiles", \
#     "/exports/nextflow/fragbuilder/analysis/molport/nodes.csv", \
#     "/exports/nextflow/fragbuilder/analysis/molport/neo"]
