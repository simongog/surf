surf
====

the SUccinct Retrival Framework.

## requirements

* gcc 4.7 or clang 4.3
* [indri](http://www.lemurproject.org/indri/) to convert indri indexes to surf input format

## installation

```
cd surf
git submodule init
git submodule update
cd surf/build
cmake ..
make
```

## building an index

```
cd surf/build
./surf_index-IDX_D -c ../collections/wikishort/
```

## querying an index

```
cd surf/build
./surf_search -c ../collections/wikishort/ -q <qryfile> -k 10
```

## creating an indri index and converting it into surf format

### create the indri index

```
cd ./indri-5.6/
./configure
make
cd buildindex
# change indri config to correct storage locations
./IndriBuildIndex ./surf/extras/gov2.indricfg
```

### convert the index into surf format

```
cd surf/tools
# change path of indri source code in Makefile
make
# ./indri_to_surf <path_to_indri_repo> <path_to_surf_repo>
./indri_to_surf ../collections/gov2indi ../collections/gov2/ 
```

## starting a daemon

an index daemon can be started in the background and listen on a specific port for search requests

```
cd build
./surf_daemon-IDX_D -c ../collections/wikishort/ -p 12345
```

## querying the search daemon

the daemon can be queried via the network or localhost

```
cd build
./surf_query -q <qry_file> -h 127.0.0.1:12345 -k 10
```

## shutting down the daemon

the daemon can be terminated via the query client.

```
cd build
./surf_query -q <qry_file> -h 127.0.0.1:12345 -k 1 -s
```














