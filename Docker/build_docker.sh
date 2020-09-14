#!/bin/bash

set -e

VERSION=`cat VERSION.txt`

docker build -t trinityctat/ctat_vif:$VERSION .
docker build -t trinityctat/ctat_vif:latest .
