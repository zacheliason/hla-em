#!/bin/bash

DEFAULT_VOLUME_PATH="/data"
VOLUME_PATH=${1:-$DEFAULT_VOLUME_PATH}
HOST_PATH=$(pwd)/src

echo "$HOST_PATH":"$VOLUME_PATH"

# Run the Docker command with the specified volume path
docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" mjinkm/work_env:latest /bin/bash
