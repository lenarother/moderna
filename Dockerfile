
FROM mattselph/ubuntu-flask

MAINTAINER Kristian Rother <krother@academis.eu>

RUN apt-get update
RUN apt-get install python-dev -f -y
RUN pip install numpy
RUN pip install biopython
RUN git clone https://github.com/lenarother/moderna.git
RUN cd moderna;python setup.py install

