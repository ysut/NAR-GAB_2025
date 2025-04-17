#!/bin/bash

set -eu

docker build -t "ps_scoring:0.1" -f psscoring/Dockerfile psscoring
docker build -t "ps_spliceai:1.3.1" -f spliceai/Dockerfile spliceai
docker build -t "ps_vep:113.4" -f vep/Dockerfile vep