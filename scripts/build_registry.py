#!/usr/bin/env python3
"""Validate every registry/clocks/*.json entry and build docs/registry.json.

Usage:
    python scripts/build_registry.py          # validate + build
    python scripts/build_registry.py --check  # validate only (used on pull requests)

Standard library only. The checks are deliberately strict about privacy: an entry
may contain model parameters and aggregate statistics, never patient-level data.
"""
import json, math, re, sys, glob, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FORBIDDEN_KEY = re.compile(r"patient|paciente|protocolo|dni|birth|nacimiento|fecha|date_of|accession|mrn", re.I)
MAX_STR = 2000         # longest string value allowed
MAX_BYTES = 100_000    # largest entry file allowed
MAX_LIST = 60          # longest list of numbers allowed (blocks pasted raw columns)
MIN_N = 100            # smallest reference / scored n accepted
STATUS = {"self-reported", "verified", "example"}
OK_FAMILY_TYPES = {"gompertz_phenoage"}


def walk(obj, path=""):
    """Yield (path, key, value) for every key in a nested JSON object."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield path, k, v
            yield from walk(v, f"{path}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, f"{path}[{i}]")


def need(d, key, kind, errs, where):
    if key not in d:
        errs.append(f"{where}: missing '{key}'")
        return None
    if not isinstance(d[key], kind):
        errs.append(f"{where}: '{key}' must be {kind}")
        return None
    return d[key]


def num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def validate(e, fname):
    errs, w = [], os.path.basename(fname)
    for path, k, v in walk(e):
        if FORBIDDEN_KEY.search(k):
            errs.append(f"{w}: key '{path}/{k}' looks like patient-level data and is not allowed")
        if isinstance(v, str) and len(v) > MAX_STR:
            errs.append(f"{w}: string at '{path}/{k}' is longer than {MAX_STR} characters")
        if isinstance(v, list) and len(v) > MAX_LIST and all(num(x) for x in v):
            errs.append(f"{w}: '{path}/{k}' is a list of {len(v)} numbers; only aggregates are allowed")
    if w != f"{e.get('id')}.json":
        errs.append(f"{w}: file name must equal id + '.json'")
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{2,60}", str(e.get("id", ""))):
        errs.append(f"{w}: id must be lowercase letters, digits and hyphens")
    if e.get("status") not in STATUS:
        errs.append(f"{w}: status must be one of {sorted(STATUS)}")
    if e.get("status") == "verified":
        errs.append(f"{w}: contributors cannot self-assign 'verified'; maintainers set it")
    lab = need(e, "lab", dict, errs, w) or {}
    need(lab, "name", str, errs, w + " lab"); need(lab, "country", str, errs, w + " lab")
    # Optional map location: a city name and its centre, never an address. At most 2 decimals (about 1 km).
    for f in ("lat", "lon"):
        if f in lab:
            v = lab[f]
            lim = 90 if f == "lat" else 180
            if not num(v) or not -lim <= v <= lim:
                errs.append(f"{w}: lab.{f} must be a number in [-{lim}, {lim}]")
            elif round(v, 2) != v:
                errs.append(f"{w}: lab.{f} may have at most 2 decimals (city centre, not an address)")
    if ("lat" in lab) != ("lon" in lab):
        errs.append(f"{w}: lab.lat and lab.lon go together")
    if "city" in lab and not isinstance(lab["city"], str):
        errs.append(f"{w}: lab.city must be a string")
    clk = need(e, "clock", dict, errs, w) or {}
    need(clk, "name", str, errs, w + " clock")
    bms = need(clk, "biomarkers", list, errs, w + " clock") or []
    model = need(clk, "model", dict, errs, w + " clock") or {}
    if model.get("type") not in OK_FAMILY_TYPES:
        errs.append(f"{w}: model.type must be one of {sorted(OK_FAMILY_TYPES)}")
    keys = set()
    for b in bms:
        for f in ("key", "label", "unit"):
            if not isinstance(b.get(f), str): errs.append(f"{w}: biomarker needs string '{f}'")
        if not num(b.get("scale")) or b["scale"] <= 0:
            errs.append(f"{w}: biomarker '{b.get('key')}' needs a finite positive 'scale'")
        if not isinstance(b.get("log"), bool): errs.append(f"{w}: biomarker '{b.get('key')}' needs boolean 'log'")
        keys.add(b.get("key"))
    betas = model.get("betas", {})
    if "age" not in betas: errs.append(f"{w}: model.betas must include 'age'")
    if set(betas) - {"age"} != keys:
        errs.append(f"{w}: model.betas keys {sorted(set(betas)-{'age'})} must match biomarker keys {sorted(keys)}")
    for f in ("rate", "m_n", "m_d", "ba_n", "ba_d", "ba_i"):
        if not num(model.get(f)): errs.append(f"{w}: model.{f} must be a number")
    for k2, v2 in betas.items():
        if not num(v2): errs.append(f"{w}: beta '{k2}' must be a number")
    rc = need(e, "recentering", dict, errs, w) or {}
    for f in ("intercept", "slope_age", "slope_male", "resid_sd"):
        if not num(rc.get(f)): errs.append(f"{w}: recentering.{f} must be a number")
    if num(rc.get("resid_sd")) and rc["resid_sd"] <= 0: errs.append(f"{w}: resid_sd must be positive")
    if not num(rc.get("reference_n")) or rc["reference_n"] < MIN_N:
        errs.append(f"{w}: recentering.reference_n must be a number >= {MIN_N}")
    cuts = rc.get("z_cutpoints", {})
    vals = [cuts.get(k) for k in ("p10", "p25", "p50", "p75", "p90")]
    if not all(num(v) for v in vals) or vals != sorted(vals):
        errs.append(f"{w}: z_cutpoints p10<p25<p50<p75<p90 must be increasing numbers")
    perf = need(e, "performance", dict, errs, w) or {}
    if not num(perf.get("n_scored")) or perf["n_scored"] < MIN_N:
        errs.append(f"{w}: performance.n_scored must be a number >= {MIN_N}")
    if not num(perf.get("r_with_age")) or not -1 <= perf["r_with_age"] <= 1:
        errs.append(f"{w}: performance.r_with_age must be in [-1, 1]")
    if not isinstance(e.get("license"), str): errs.append(f"{w}: missing license for the entry")
    return errs


def main():
    check_only = "--check" in sys.argv
    entries, errors = [], []
    for f in sorted(glob.glob(os.path.join(ROOT, "registry", "clocks", "*.json"))):
        try:
            if os.path.getsize(f) > MAX_BYTES:
                raise ValueError(f"file larger than {MAX_BYTES} bytes")
            def _no_const(c):
                raise ValueError(f"'{c}' is not valid JSON")
            e = json.load(open(f, encoding="utf-8"), parse_constant=_no_const)
        except Exception as ex:
            errors.append(f"{os.path.basename(f)}: invalid JSON ({ex})"); continue
        errors += validate(e, f)
        entries.append(e)
    ids = [e.get("id") for e in entries]
    if len(ids) != len(set(ids)): errors.append("duplicate ids in registry")
    if errors:
        print("Registry validation FAILED:"); [print("  -", x) for x in errors]; sys.exit(1)
    print(f"Registry OK: {len(entries)} entries")
    if not check_only:
        out = os.path.join(ROOT, "docs", "registry.json")
        json.dump({"schema": 1, "clocks": entries}, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        print("wrote", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
