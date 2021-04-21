#!/bin/bash

# install twolevel_ants_dbm
git clone https://github.com/CoBrALab/twolevel_ants_dbm.git $PWD/twolevel_ants_dbm
cd $PWD/twolevel_ants_dbm
git checkout 3623b5bf683473887b7d63c1a1b24c0362ed8396

# install twolevel_ants_dbm
git clone https://github.com/CoBrALab/minc-toolkit-extras.git $PWD/minc-toolkit-extras
cd $PWD/minc-toolkit-extras
git checkout 3dec9b81b0e59c7daa90ea1901469111b2374182
