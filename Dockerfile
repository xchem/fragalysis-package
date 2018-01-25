FROM informaticsmatters/rdkit:Release_2017_03_3
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
RUN git clone https://github.com/rdkit/mmpdb /usr/local/mmpdb
RUN pip install /usr/local/mmpdb 
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
