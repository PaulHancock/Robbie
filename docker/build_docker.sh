#! /usr/bin/env bash
cd ..
version=$(robbie_version.sh)

# build the container and update tag
docker build . -f docker/Dockerfile -t "paulhancock/robbie-next:${version}" && \
docker tag paulhancock/robbie-next:${version} paulhancock/robbie-next:latest
cd -