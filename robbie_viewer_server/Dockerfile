FROM python:3.9

LABEL maintainer="Paul Hancock <paul.hancock@curtin.edu.au>"

# Set the Working directory
WORKDIR /robbie_viewer

# Copy the viewer
COPY . .

# Run pip install for requirements
RUN pip install -r requirements.txt

# EXPOSE PORTS
EXPOSE 5006

CMD bokeh serve . --dev --args ./results