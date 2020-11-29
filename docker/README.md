## Build image

```bash
docker build -f Dockerfile -t kmtricks-d .
```

## Run
Entrypoint corresponds to [kmtricks.py](kmtricks.py) pipeline, so you can use:

```bash
docker run --rm -i -t -v $PWD/SHARED:/tmp kmtricks-d [kmtricks cli args]
```
**Example:**

```bash
docker run --rm -i -t -v $PWD/SHARED:/tmp kmtricks-d run --run-dir /tmp/RUN_DIR --file /tmp/data/fof.txt 
```

```
./SHARED/
├── data
│   ├── 1.fasta
│   ├── 2.fasta
│   └── fof.txt
└── RUN_DIR        // create by kmtricks
```
