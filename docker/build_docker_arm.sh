#! /usr/bin/env bash
cd ..
version=$(robbie_version.sh)

# build the container and update tag
docker buildx build . --platform linux/amd64 --no-cache -f docker/Dockerfile -t "paulhancock/robbie-next:${version}" && \
docker tag paulhancock/robbie-next:${version} paulhancock/robbie-next:latest
cd -