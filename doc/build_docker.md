# Building megSAP with Docker

## Instructions
The build instructions are quite profound:

1. Install Docker
2. Start it!
3. Check out the repo (`git clone https://github.com/imgag/megSAP.git`)
4. `cd megSAP`
5. `docker build -t "megSAP:local" .`

### Building for Ubuntu 18.04
The docker image supports for building with different platforms, particularily Ubuntu 16 / 18.

To do so, supply the `UBUNTU_VERSION` build argument (either 16 or 18)

```
docker build --build-arg UBUNTU_VERSION=18 -t megSAP:local .
```

### Uploading Docker to the registry
- [ ] TODO
