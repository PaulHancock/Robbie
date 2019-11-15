FROM python:3.6

MAINTAINER Paul Hancock


RUN apt update
# install aegean tools
RUN pip install aegeantools
RUN pip install pandas

# install fits_warp
RUN pip install psutil
RUN mkdir /tmp/build && cd /tmp/build && git clone https://github.com/nhurleywalker/fits_warp.git && cd fits_warp\
    && git checkout py3 && chmod ugo+x fits_warp.py && mv fits_warp.py /usr/bin/.

# 'install' the stilts/topcat file
RUN cd /usr/local/lib && wget http://www.star.bris.ac.uk/~mbt/stilts/stilts.jar

RUN apt install -y openjdk-11-jdk