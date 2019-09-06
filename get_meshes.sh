#!/bin/bash

curl -L https://www.dropbox.com/s/x6aqbub3qjqd527/meshes.zip?dl=0 > meshes.zip
unzip meshes.zip
rm meshes.zip
rm -r -f __MACOSX/
