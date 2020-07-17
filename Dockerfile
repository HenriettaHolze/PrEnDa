FROM continuumio/miniconda3:latest

COPY . /app

RUN conda env create -f /app/env.yaml
RUN source activate PrEnDa

RUN python /app/dashboard_local.py

CMD ["/bin/bash"]