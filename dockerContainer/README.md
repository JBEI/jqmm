# Installing all jQMM dependencies using Docker 
Docker provides a way to **increase reproducibility** for our jQMM results. Docker containers (https://www.docker.com/) wrap a piece of sofware in a complete file system containing everything needed to run: code, systems libraries, systems tools and anything else that can be installed in a server. This guarantees that it will always run correctly and in the same way, regardless of the system environment it is running in. The jQMM docker container can be run on any server or cloud computing service such as Amazon Web Services (AWS), Google Cloud Platform or Microsoft Azure. This container does **NOT** include GAMS, CPLEX or CONOPT licenses, which must be provided by the user.

Here are the instructions for a jQMM installation through docker:
## 0) Install docker:

https://docs.docker.com/engine/installation/

## 1) Launch container:

docker run --name jqmm -p 8888:8888 -d mhgarci1/jqmm

## 2) Install GAMS:

docker exec -it jqmm ./install_gams.sh

## 3) Copy in the gamslice.txt license from the current folder.
***NOTICE:*** This command will only work if you have a **valid GAMS license**, and run the command from within a folder already containing that license, which must be named "gamslice.txt". If you do not have a valid GAMS license, this step can be skipped and GAMS will operate in trial mode, where only small models can be solved.

docker cp gamslice.txt jqmm:/gams/gamslice.txt

## 4) Run the GAMS script which finalizes the GAMS install. 
**Ignore the instructions and questions on the screen**, as our workflow performs those steps automatically.

docker exec -it jqmm bash -c 'cd /gams; yes 1 | ./gamsinst'

## 5) Stop container, restart interactively, and  manually open a browser to the url it suggests:

docker stop jqmm

docker start -i jqmm

## Optional
### To **completely uninstall jQMM** run the following three commands.

docker stop jqmm

docker rm jqmm

docker rmi mhgarci1/jqmm

### To build the docker image locally (only for mhgarci1 docker account owner)

docker build -t mhgarci1/jqmm:latest .

docker login

docker push mhgarci1/jqmm


