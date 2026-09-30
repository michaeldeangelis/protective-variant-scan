"""Normalized burden-table schema, prereg config loading, table validation."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

COLUMNS = ["gene", "trait", "mask", "cohort", "beta", "se", "p", "n_carriers", "n_total", "source"]
KEY = ["gene", "trait", "mask", "cohort"]
MASKS = ("plof", "dmis", "syn")
COHORTS = ("discovery", "replication", "discovery_eur")
DOMAINS = ("cognitive", "physical", "metabolic")


# C10c: sha256 of the committed config/prereg.yaml. Changes only with a dated ledger entry.
PINNED_CONFIG_SHA256 = "be69827539daac4c78b4c8f6ac7d931723abb41ab0e739c0505c357689310a9f"


@dataclass(frozen=True)
class Config:
    sha256: str
    ledger_entry: str
    panel_sign: dict
    panel_domain: dict
    tradeoff_sign: dict
    control_sign: dict
    discovery_p: float
    replication_p: float
    proxies: dict          # panel trait -> {"trait": str, "benefit": int}
    independent_sources: tuple  # source-name tokens accepted as independent of UK Biobank
    ukb_sources: tuple     # source-name tokens that mark a row as UK Biobank (overrides the allow list)
    tradeoff_p: float
    min_tradeoffs_screened: int
    lambda_gc_max: float
    syn_contaminated_max_fraction: float
    syn_min_coverage: float
    syn_min_rows: int
    syn_max_p1_fraction: float
    replication_sign_max_p: float
    replication_sign: tuple
    positive_controls: tuple
    qualifying_domains: tuple
    lipid_genes: frozenset

    def sign(self, trait: str) -> int:
        """Direction of benefit (+1/-1) for any declared trait, proxy included."""
        for d in (self.panel_sign, self.tradeoff_sign, self.control_sign):
            if trait in d:
                return d[trait]
        for px in self.proxies.values():
            if px["trait"] == trait:
                return px["benefit"]
        raise KeyError(f"no declared direction of benefit for trait {trait!r}")

    @property
    def panel(self) -> tuple:
        return tuple(self.panel_sign)

    @property
    def tradeoff(self) -> tuple:
        return tuple(self.tradeoff_sign)


def _signs(d: dict, where: str) -> dict:
    out = {}
    for k, v in d.items():
        if v["benefit"] not in (1, -1):
            raise ValueError(f"{where}.{k}.benefit must be +1 or -1")
        out[k] = int(v["benefit"])
    return out


def load_config(path) -> Config:
    raw_bytes = Path(path).read_bytes()
    y = yaml.safe_load(raw_bytes)
    panel_sign = _signs(y["panel"], "panel")
    panel_domain = {k: v["domain"] for k, v in y["panel"].items()}
    tradeoff_sign = _signs(y["tradeoff"], "tradeoff")
    control_sign = _signs(y["control_traits"], "control_traits")
    proxies = {k: {"trait": v["trait"], "benefit": int(v["benefit"])} for k, v in y["replication"]["proxies"].items()}

    d = y["discovery"]
    if len(panel_sign) != d["n_traits"]:
        raise ValueError(f"panel has {len(panel_sign)} traits, discovery.n_traits={d['n_traits']}")
    formula = d["alpha"] / (d["n_genes"] * d["n_traits"])
    if abs(d["p_threshold"] - formula) / formula > 0.02:
        raise ValueError(f"discovery.p_threshold {d['p_threshold']} inconsistent with alpha/(n_genes*n_traits)={formula}")
    if any(v not in DOMAINS for v in panel_domain.values()):
        raise ValueError("panel domain must be one of " + ", ".join(DOMAINS))
    if not set(proxies) <= set(panel_sign):
        raise ValueError("replication.proxies keys must be panel traits")
    if any(px["benefit"] not in (1, -1) for px in proxies.values()):
        raise ValueError("proxy benefit must be +1 or -1")
    known = set(panel_sign) | set(tradeoff_sign) | set(control_sign)
    pos = tuple(y["controls"]["positive"])
    for c in pos:
        if c["trait"] not in known:
            raise ValueError(f"positive control {c['id']}: unknown trait {c['trait']}")
    min_screened = int(y["tradeoff_screen"]["min_tradeoffs_screened"])
    if not 1 <= min_screened <= len(tradeoff_sign):
        raise ValueError("tradeoff_screen.min_tradeoffs_screened must be between 1 and N_tradeoff")
    if not 0 < float(y["controls"]["syn_contaminated_max_fraction"]) < 1:
        raise ValueError("controls.syn_contaminated_max_fraction must be in (0, 1)")
    for c in y["controls"]["replication_sign"]:
        if c["expected_beta_sign"] not in (1, -1):
            raise ValueError(f"replication_sign {c['id']}: expected_beta_sign must be +1 or -1")
    q = tuple(y["verdict"]["qualifying_domains"])
    if not set(q) <= set(DOMAINS):
        raise ValueError("verdict.qualifying_domains invalid")

    return Config(
        sha256=hashlib.sha256(raw_bytes).hexdigest(),
        ledger_entry=y["ledger_entry"],
        panel_sign=panel_sign,
        panel_domain=panel_domain,
        tradeoff_sign=tradeoff_sign,
        control_sign=control_sign,
        discovery_p=float(d["p_threshold"]),
        replication_p=float(y["replication"]["one_sided_p"]),
        proxies=proxies,
        independent_sources=tuple(t.lower() for t in y["replication"]["independent_sources"]),
        ukb_sources=tuple(t.lower() for t in y["replication"]["ukb_overlapping_sources"]),
        tradeoff_p=float(y["tradeoff_screen"]["alpha"]) / len(tradeoff_sign),
        min_tradeoffs_screened=min_screened,
        lambda_gc_max=float(y["controls"]["lambda_gc_max"]),
        syn_contaminated_max_fraction=float(y["controls"]["syn_contaminated_max_fraction"]),
        syn_min_coverage=float(y["controls"]["syn_min_coverage"]),
        syn_min_rows=int(y["controls"]["syn_min_rows"]),
        syn_max_p1_fraction=float(y["controls"]["syn_max_p1_fraction"]),
        replication_sign_max_p=float(y["controls"]["replication_sign_max_p"]),
        replication_sign=tuple(y["controls"]["replication_sign"]),
        positive_controls=pos,
        qualifying_domains=q,
        lipid_genes=frozenset(g.upper() for g in y["verdict"]["lipid_pathway_genes"]),
    )


def validate(df: pd.DataFrame) -> pd.DataFrame:
    """Return a cleaned copy of a normalized table; raise ValueError on schema violations."""
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"missing columns: {missing}")
    df = df[COLUMNS].copy()
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    df["trait"] = df["trait"].astype(str).str.strip().str.lower()
    df["mask"] = df["mask"].astype(str).str.strip().str.lower()
    df["cohort"] = df["cohort"].astype(str).str.strip().str.lower()
    for c in ("beta", "se", "p", "n_carriers", "n_total"):
        df[c] = pd.to_numeric(df[c], errors="raise")
    bad = sorted(set(df["mask"]) - set(MASKS))
    if bad:
        raise ValueError(f"unknown mask values: {bad}")
    bad = sorted(set(df["cohort"]) - set(COHORTS))
    if bad:
        raise ValueError(f"unknown cohort values: {bad}")
    n0 = len(df)
    df = df.dropna(subset=["beta", "p"])          # se may be undefined (beta = 0 or p = 1 rows, C10d)
    dropped = n0 - len(df)
    if ((df["p"] < 0) | (df["p"] > 1)).any():
        raise ValueError("p outside [0, 1]")
    if (df["se"] < 0).any():
        raise ValueError("negative se")
    dup = df.duplicated(KEY, keep=False)
    if dup.any():
        ex = df.loc[dup, KEY].head(3).to_dict("records")
        raise ValueError(f"{int(dup.sum())} rows share a (gene, trait, mask, cohort) key, e.g. {ex}")
    df = df.reset_index(drop=True)
    df.attrs["dropped_nan_rows"] = dropped
    return df


def load_burden(data_dir) -> pd.DataFrame:
    """Concatenate and validate every top-level data_dir/burden*.csv[.gz]."""
    files = sorted(Path(data_dir).glob("burden*.csv")) + sorted(Path(data_dir).glob("burden*.csv.gz"))
    if not files:
        raise FileNotFoundError(f"no burden*.csv[.gz] files in {data_dir}")
    df = validate(pd.concat([pd.read_csv(f, dtype={"gene": str, "trait": str}) for f in files], ignore_index=True))
    synth = {x for x in df["source"].astype(str).unique() if x.lower().startswith("synthetic")}
    real = set(df["source"].astype(str).unique()) - synth
    if synth and real:
        raise ValueError(f"{data_dir} mixes synthetic sources {sorted(synth)} with non-synthetic {sorted(real)}; refusing to run")
    df.attrs["files"] = [f.name for f in files]
    return df
