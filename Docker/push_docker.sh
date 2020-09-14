#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker push trinityctat/ctat_vif:$VERSION
docker push trinityctat/ctat_vif:latest
