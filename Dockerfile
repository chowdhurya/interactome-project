FROM continuumio/miniconda

ADD . /code
WORKDIR /code

RUN conda install --file conda-requirements.txt -y
