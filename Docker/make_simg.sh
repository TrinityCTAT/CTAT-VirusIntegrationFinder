s#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build -F ctat_vif.v${VERSION}.simg docker://trinityctat/ctat_vif:$VERSION

singularity exec -e ctat_vif.v${VERSION}.simg ctat-VIF.py --help


