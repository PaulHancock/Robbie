FROM python:3.6

LABEL maintainer="Paul Hancock <paul.hancock@curtin.edu.au>"

# non-python dependencies
RUN apt update && \
    apt install -y openjdk-11-jdk sqlite3 && \
    apt-get autoremove -y && \
    apt-get clean

# stilts 
RUN cd /usr/local/lib && wget http://www.star.bris.ac.uk/~mbt/stilts/stilts.jar && \
    cd /usr/local/bin && echo 'java -jar /usr/local/lib/stilts.jar "$@"' > /usr/local/bin/stilts && chmod ugo+x /usr/local/bin/stilts

# python modules
RUN pip install aegeantools pandas psutil && \
    mkdir /tmp/build && cd /tmp/build && \
    git clone https://github.com/nhurleywalker/fits_warp.git && cd fits_warp && \
    git checkout py3 && chmod ugo+x fits_warp.py && mv fits_warp.py /usr/bin/. 

# add Robbie scripts (always)
ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN cd /tmp/build && git clone https://github.com/PaulHancock/Robbie.git && \
    cd Robbie && chmod ugo+x *.py && mv *.py /usr/bin/. && \
    rm -rf /tmp/build
