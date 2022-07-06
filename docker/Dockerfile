FROM python:3.9

LABEL maintainer="Paul Hancock <paul.hancock@curtin.edu.au>"

# non-python dependencies
RUN apt update && \
    apt install -y openjdk-11-jdk swarp && \
    apt-get autoremove -y && \
    apt-get clean

# stilts
RUN cd /usr/local/lib && wget http://www.star.bris.ac.uk/~mbt/stilts/stilts.jar && \
    cd /usr/local/bin && echo 'java -jar /usr/local/lib/stilts.jar "$@"' > /usr/local/bin/stilts && chmod ugo+x /usr/local/bin/stilts

# python dependencies, with specific versions specified for longevity, and Robbie scripts (always)
WORKDIR /tmp/build
ADD . /tmp/build
# using pip install . will break the shebang lines of some scripts so stick with python setup.py install
RUN pip install -r requirements.txt && \
    python setup.py install && \
    rm -rf /tmp/build

#  set env
ENV HOME=/tmp