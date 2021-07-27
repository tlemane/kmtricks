# Build image

```bash
git clone https://github.com/tlemane/kmtricks
cd kmtricks/docker
docker build -f Dockerfile -t kmtricks-d
```

# Run

```bash
docker run --rm -i -t -v $PWD/SHARED:/tmp kmtricks-d <kmtricks cli args>
```

The default entrypoint corresponds to `kmtricks`. To run `kmtricks-socks`, use:

```bash
docker run --rm -i -t -v $PWD/SHARED:/tmp --entrypoint /opt/kmtricks/bin/kmtricks-socks kmdiff-d <kmtricks cli args>
```
