## queries

This repository contains the SPARQL queries (one per file) to generate
a [grlc](https://github.com/CLARIAH/grlc)-based Web API that exposes data in the
[pbg-ld](https://github.com/candYgene/pbg-ld) platform.


**Web API endpoint(s)**
- base URL `http://CNAME:PORT` (default: `http://localhost:8088) followed by
- remote path `/api-git/candYgene/queries/` or
  - requires `GRLC_GITHUB_ACCESS_TOKEN` to be set
- local path `/api-local/`
  - requires `docker cp queries grlc:/home/grlc/`
