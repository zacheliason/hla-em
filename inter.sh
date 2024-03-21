#!/bin/bash

VOLUME_PATH="/data"
HOST_PATH=$(pwd)/hla-em

echo "$HOST_PATH":"$VOLUME_PATH"
echo ""

# Run the Docker command with the specified volume path
docker run -it --rm -v "$HOST_PATH":"$VOLUME_PATH" mjinkm/work_env:latest /bin/bash
