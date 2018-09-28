# Installing and running via Docker

If you haven't done so already, please have a look at [Docker](https://docker.io) and [install it for your platform](https://docs.docker.com/install/).

After you've successfully installed Docker you are ready to run megSAP.

## Running megSAP with Docker

The most basic usage of _megSAP_ via Docker is

```
docker run -it imgag/megSAP ls src/Pipelines
```

This will download the _megSAP_ container and list the available pipelines. Please refer to the [documentation](https://github.com/imgag/megSAP#documentation) on how to use individual pipelines.

***Note:*** Please bear in mind that megSAP uses a very versatile set of tools. Because of this the container is rather large and the initial download may be 1GB or larger, thus taking some time.

### Using a custom configuration file

If you would like to pass a custom configuration file (_settings.ini_) you can do it like so:

```
docker run -it -v /path/to/your/settings.ini:/home/ubuntu/megSAP/settings.ini imgag/megSAP ls pipelines
```

Generally the `data` and `src` folders of megSAP are available with all components installed, so you can [mount those](https://docs.docker.com/storage/volumes/) at any time.

### Running megSAP on a more recent Ubuntu
It is possible to build megSAP on a recent Ubuntu. Please refer to the [docker build guide](/build_docker.md).
