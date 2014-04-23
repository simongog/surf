surf
====

# installation

```
cd surf/build
cmake ..
make
```

# building an index

```
cd surf/build
./surf_index -c ../collections/wikishort/
```

# querying an index

```
cd surf/build
./surf_search -c ../collections/wikishort/ -q <qryfile> -k 10
```

# creating an indri index and converting it into surf format

create the indri index

```
cd ./indri-5.6/
./configure
make
cd buildindex
./IndriBuildIndex ./surf/extras/gov2.indriconfig
```

convert the index into surf format

```
cd surf/tools
# change path of indri source code in Makefile
make
# ./indri_to_surf <path_to_indri_repo> <path_to_surf_repo>
./indri_to_surf ../collections/gov2indi ../collections/gov2/ 
```



