#!/bin/bash

DEFAULT_VOLUME_PATH="/data"
VOLUME_PATH=${1:-$DEFAULT_VOLUME_PATH}
HOST_PATH=$(pwd)/src

echo "mounting $HOST_PATH":"$VOLUME_PATH"

# Run the Docker command with the specified volume path
docker run --rm -v "$HOST_PATH":"$VOLUME_PATH" hla_em:latest /bin/bash -c "$VOLUME_PATH"/run_hla.sh $VOLUME_PATH
