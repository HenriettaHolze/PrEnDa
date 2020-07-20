FROM continuumio/miniconda3:latest
RUN conda update conda

RUN mkdir /app
COPY env.yaml /app/

RUN conda env create -f /app/env.yaml

COPY . /app

RUN apt-get -y -q update
RUN apt-get -y -q install curl
RUN apt-get -y -q install libsm6 libxrender1 libfontconfig1

CMD ["/bin/bash", "-c", ". /app/init.sh; python /app/dashboard_local.py"]