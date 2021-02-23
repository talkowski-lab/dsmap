# Building DSMap Docker images  

All Docker images used in this project can be built with `build_docker.py`.  

Docker images included in this script are:  
  * `us.gcr.io/broad-dsmap/athena`
  * `us.gcr.io/broad-dsmap/athena-cloud`
  * `us.gcr.io/broad-dsmap/dsmap`

Example usage to rebuild all Docker images:
```
./build_docker.py \
	--images All \
	--athena-hash "$desired_athena_github_commit_hash" \
	--dsmap-hash "$desired_dsmap_github_commit_hash" \
	--tag "my_new_image" \
	--gcr-project broad-dsmap \
	--update-latest
```  
