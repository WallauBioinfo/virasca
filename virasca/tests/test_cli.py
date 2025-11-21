import os
import subprocess
import sys
from pathlib import Path

def test_cli_help():
    print("Testing CLI help...")
    result = subprocess.run(["python3", "virasca/virasca/cli.py", "--help"], capture_output=True, text=True)
    if result.returncode != 0:
        print("CLI help failed")
        print(result.stderr)
        sys.exit(1)
    print("CLI help passed")

def test_configure_database_dryrun():
    print("Testing configure-database dry-run...")
    # We need to ensure we are in the right directory or set python path
    env = os.environ.copy()
    env["PYTHONPATH"] = os.getcwd() + "/virasca"
    
    # Running from repo root, but setting cwd to virasca/virasca (package root?)
    # Actually, let's run from repo root and point to cli.py correctly.
    # Repo root: /Users/filipedezordi/Git/virasca
    # cli.py: virasca/virasca/cli.py
    
    result = subprocess.run(
        ["python3", "virasca/virasca/cli.py", "configure-database", "--dry-run"],
        capture_output=True,
        text=True,
        env=env,
        cwd=os.getcwd() # Run from repo root
    )
    if result.returncode != 0:
        print("configure-database dry-run failed")
        print(result.stderr)
    else:
        print("configure-database dry-run passed")
        print(result.stdout)

def test_run_dryrun():
    print("Testing run dry-run...")
    env = os.environ.copy()
    env["PYTHONPATH"] = os.getcwd() + "/virasca"
    
    # Create dummy files
    Path("test.fa").touch()
    Path("r1.fq").touch()
    Path("r2.fq").touch()
    
    try:
        result = subprocess.run(
            ["python3", "virasca/virasca/cli.py", "run", "--input", "test.fa", "--reads-r1", "r1.fq", "--reads-r2", "r2.fq", "--dry-run"],
            capture_output=True,
            text=True,
            env=env,
            cwd=os.getcwd()
        )
        if result.returncode != 0:
            print("run dry-run failed")
            print(result.stderr)
        else:
            print("run dry-run passed")
            print(result.stdout)
    finally:
        # Cleanup
        if os.path.exists("test.fa"): os.remove("test.fa")
        if os.path.exists("r1.fq"): os.remove("r1.fq")
        if os.path.exists("r2.fq"): os.remove("r2.fq")

if __name__ == "__main__":
    test_cli_help()
    # We need to be careful about CWD.
    # The cli.py is in virasca/virasca/cli.py
    # The Snakefile is in virasca/virasca/workflow/Snakefile
    # The config is in virasca/virasca/config/config.yaml
    # If we run from repo root:
    # python3 virasca/virasca/cli.py
    # Inside cli.py: get_snakefile returns .../workflow/Snakefile
    # Snakemake runs with that snakefile.
    # Snakemake resolves relative paths in Snakefile relative to Snakefile? Or CWD?
    # "configfile: config/config.yaml" -> relative to CWD usually.
    # So we should run from virasca/virasca directory?
    # Or we should make paths absolute in Snakefile.
    # Let's try running the tests and see.
    test_configure_database_dryrun()
    test_run_dryrun()
