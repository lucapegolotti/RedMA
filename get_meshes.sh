#!/bin/bash

# curl -L https://www.dropbox.com/s/05jy03munytgwpv/meshes.zip?dl=0 > meshes.zip
curl -L https://www.dropbox.com/sh/aqtk18jekylksya/AADQlM74JjuGEhY-AHizUo27a?dl=0 > meshes.zip
mkdir meshes
mv meshes.zip meshes
cd meshes
unzip meshes.zip
rm meshes.zip
rm -r -f __MACOSX/
