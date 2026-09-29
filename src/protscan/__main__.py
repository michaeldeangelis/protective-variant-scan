import argparse
import sys

from .run import run_pipeline


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="protscan")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the preregistered scan")
    r.add_argument("--config", required=True)
    r.add_argument("--data", required=True)
    r.add_argument("--out", required=True, help="results JSON path; report.md is written beside it")
    a = ap.parse_args(argv)
    res = run_pipeline(a.config, a.data, a.out)
    v = res["verdict"]
    print(f"verdict: {v['label']} ({v['reason']})")
    print(f"controls valid: {res['controls']['valid']}")
    print(f"wrote {a.out} and report.md")
    if not res["config_matches_ledger"]:
        print("NON-PREREGISTERED: config sha256 does not match the pinned ledger value", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
