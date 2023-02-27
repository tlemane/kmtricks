# Build image

```bash
git clone https://github.com/tlemane/kmtricks
cd kmtricks/docker
docker build -f Dockerfile -t kmtricks-d .
```

# Run

```bash
docker run --rm -v $PWD/SHARED:/tmp kmtricks-d <kmtricks cli args>
```

The default entrypoint corresponds to `kmtricks`. To run `kmtricks` with plugin support, use:

```bash
docker run --rm -v $PWD/SHARED:/tmp --entrypoint kmtricksp kmtricks-d <kmtricks cli args>
```
