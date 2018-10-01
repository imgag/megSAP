# Installing and running via Docker

If you haven't done so already, please have a look at [Docker](https://docker.io) and [install it for your platform](https://docs.docker.com/install/).

After you've successfully installed Docker you are ready to run megSAP.

## Running megSAP with Docker

To set up megSAP using docker execute the below commands

```
docker run -it -v ./genomes:/megSAP/data/genomes -v imgag/megSAP cd data && ./download_GRCh37.sh
docker run -it -v ./genomes:/megSAP/data/genomes -v ./dbs:/megSAP/data/dbs imgag/megSAP cd data && ./download_dbs
```

This will download and execute the _megSAP_ container, thus downloading the needed databases in your **current working directory**.
After setup you can use the _megSAP_ container to run pipelines. Please refer to the [documentation](./pipelines.md) on how to use individual pipelines.

***Note:*** Please bear in mind that _megSAP_ uses a very versatile set of tools. Because of this the container is rather large and the initial download may be 1GB or larger, thus taking some time.

### Using a custom configuration file

If you would like to pass a custom configuration file (_settings.ini_) you can do it like so:

```
docker run -it -v /path/to/your/settings.ini:/home/ubuntu/megSAP/settings.ini imgag/megSAP ls pipelines
```

Generally the `data` and `src` folders of megSAP are available with all components installed, so you can [mount those](https://docs.docker.com/storage/volumes/) at any time.

### Running megSAP on a more recent Ubuntu
It is possible to build megSAP on a recent Ubuntu. Please refer to the [docker build guide](./build_docker.md).
