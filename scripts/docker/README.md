# Building DSMap Docker images  

All Docker images used in this project can be built with `build_docker.py`.  

Docker images included in this script are:  
  * `us.gcr.io/broad-dsmap/athena`
  * `us.gcr.io/broad-dsmap/athena-cloud`
  * `us.gcr.io/broad-dsmap/dsmap`
  * `us.gcr.io/broad-dsmap/dsmap-cromwell`
  * `us.gcr.io/broad-dsmap/dsmap-r`

Example usage to rebuild all Docker images using most recent repo commits:
```
./build_docker.py \
	--images all \
	--tag "my_new_image" \
	--gcr-project broad-dsmap \
	--update-latest
```  
