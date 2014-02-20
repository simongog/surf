surf
====

SUccinct Retrieval Framework

## Directory structure

    * build: 
    * collections: contains collection data;
        each collection in its own subdirectory.
    * external: Contains external libraries
    * include/surf: Contains headers.
    * src: Contains surf sources.
      - surf_index.cpp
      - surf_search.cpp
      - surf_profile.cpp
    * tools: Contains source of tools. E.g.
      - qry2intqry
      - indri2surf converter

## Installation

    * git clone git@github.com:simongog/surf.git
    * cd surf
    * git submodule init
    * git submodule update --recursive
    * cd build
    * cmake ..

