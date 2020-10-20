#!/bin/bash

curl -L https://www.dropbox.com/s/05jy03munytgwpv/meshes.zip?dl=0 > meshes.zip
unzip meshes.zip
rm meshes.zip
rm -r -f __MACOSX/
