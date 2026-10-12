import gzip

import pytest


@pytest.fixture
def write_gz(tmp_path):
    """write_gz(name, header, rows, preamble=()) -> path of a gzipped TSV."""

    def _write(name, header, rows, preamble=()):
        path = tmp_path / name
        lines = list(preamble) + ["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]
        with gzip.open(path, "wt") as fh:
            fh.write("\n".join(lines) + "\n")
        return path

    return _write
