"""Minimum free energy and secondary structure of the guide RNA.

Two backend, same return shape:

* local, when the ViennaRNA Python binding is importable (`import RNA`);
* the RNAfold web server, otherwise.

The web server takes exactly one sequence per job: a multi record FASTA is
silently concatenated into a single sequence, which would return a quietly
wrong answer, so one job is submitted per distinct guide. A JSON cache on
disk keyed by sequence, temperature and scaffold makes a repeat run free.
"""

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

#: FASTA header sent to the server. The result file is named after it.
_JOB_NAME = "guide"
_JOB_PATTERN = re.compile(r"/RNAfold/([A-Za-z0-9_-]+)")
_RESULT_PATTERN = re.compile(r"^([.()\[\]{}<>]+)\s+\(\s*(-?\d+\.\d+)\)\s*$", re.M)
_AGENT = {"User-Agent": "gRNADesign (iGEM, academic use)"}

FAILED = (None, float("nan"))


def local_available():
    try:
        import RNA  # noqa: F401
    except ImportError:
        return False
    return True


def _fold_local(sequence, temperature):
    import RNA

    model = RNA.md()
    model.temperature = float(temperature)
    fold_compound = RNA.fold_compound(sequence, model)
    structure, energy = fold_compound.mfe()
    return structure, float(energy)


def _fold_server(sequence, url, temperature, timeout=240.0, poll=2.0):
    """Submit one sequence, then poll the job's plain text result file.

    The result file 404s while the job is queued and becomes readable the
    moment it is done, which is a cleaner completion signal than scraping the
    status page.
    """
    payload = urllib.parse.urlencode({
        "PAGE": "2",
        "SCREEN": f">{_JOB_NAME}\n{sequence}",
        "method": "mfe",
        "noLP": "on",
        "dangling": "d2",
        "param": "rna2004",
        "Temp": f"{float(temperature):g}",
        "salt": "1.021",
        "proceed": "",
    }).encode()

    request = urllib.request.Request(url, data=payload, headers=_AGENT)
    with urllib.request.urlopen(request, timeout=60) as response:
        page = response.read().decode("utf-8", "replace")

    found = _JOB_PATTERN.search(page)
    if not found:
        raise RuntimeError("the RNAfold server did not return a job identifier")
    result_url = urllib.parse.urljoin(
        url, f"/RNAfold/{found.group(1)}/{_JOB_NAME}.out")

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(
                urllib.request.Request(result_url, headers=_AGENT), timeout=60
            ) as response:
                text = response.read().decode("utf-8", "replace")
            break
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
            time.sleep(poll)
    else:
        raise TimeoutError(f"the RNAfold job did not finish within {timeout:g}s")

    match = _RESULT_PATTERN.search(text)
    if not match:
        raise RuntimeError(f"could not read an MFE out of {text!r}")
    return match.group(1), float(match.group(2))


def fold(sequence, *, backend="auto", url, temperature=37.0, worker=4,
         cache_path=None, scaffold="", log=print):
    """Fold every sequence, returning {sequence: (structure, mfe)}.

    `sequence` is the bare spacer. `scaffold` is appended before folding, so
    the structure reported covers whatever the caller asked to fold. A guide
    that could not be folded maps to (None, nan) instead of failing the run.
    """
    distinct = sorted(set(sequence))
    if not distinct:
        return {}

    if backend == "auto":
        backend = "local" if local_available() else "server"
    if backend == "local" and not local_available():
        raise RuntimeError(
            "mfe_backend is 'local' but the ViennaRNA Python binding is not "
            "importable; install it (pip install ViennaRNA) or use 'server'"
        )

    cache_path = Path(cache_path) if cache_path else None
    cache = {}
    if cache_path and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except (OSError, ValueError):
            cache = {}

    def key(spacer):
        return f"{spacer}|{scaffold}|{float(temperature):g}"

    pending = [spacer for spacer in distinct if key(spacer) not in cache]
    log(f"  {len(distinct)} distinct guide, {len(distinct) - len(pending)} cached, "
        f"{len(pending)} to fold via {backend}")

    def run(spacer):
        target = (spacer + scaffold).replace("T", "U")
        if backend == "local":
            return spacer, _fold_local(target, temperature)
        return spacer, _fold_server(target, url, temperature)

    failure = 0
    if pending:
        worker = 1 if backend == "local" else max(1, int(worker))
        with ThreadPoolExecutor(max_workers=worker) as pool:
            for done, (spacer, result) in enumerate(
                _iter_result(pool, run, pending), start=1
            ):
                if result is None:
                    failure += 1
                else:
                    cache[key(spacer)] = list(result)
                if done % 10 == 0 or done == len(pending):
                    log(f"  folded {done}/{len(pending)}")

    if cache_path and pending:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=0, sort_keys=True))

    if failure:
        log(f"  {failure} guide could not be folded; its MFE is left empty")

    return {
        spacer: tuple(cache[key(spacer)]) if key(spacer) in cache else FAILED
        for spacer in distinct
    }


def _iter_result(pool, run, pending):
    """Yield (sequence, result) in completion order, turning a per sequence
    failure into a None result rather than losing the whole batch."""
    future = {pool.submit(run, spacer): spacer for spacer in pending}
    for task in as_completed(future):
        try:
            yield task.result()
        except Exception:  # noqa: BLE001 - one bad fold must not end the run
            yield future[task], None
