#! /usr/bin/env python
"""
Setup for Robbie
"""
import os
from setuptools import setup
from subprocess import check_output


# The following two functions were taken from the repo: https://github.com/pyfidelity/setuptools-git-version/blob/master/setuptools_git_version.py

def format_version(version, fmt="{tag}_{gitsha}"):
    parts = version.split("-")
    if len(parts) < 4:
        return parts[0]
    assert len(parts) in (3, 4)
    dirty = len(parts) == 4
    tag, count, sha = parts[:3]
    if count == "0" and not dirty:
        return tag
    return fmt.format(tag=tag, commitcount=count, gitsha=sha.lstrip("g"))


def get_git_version():
    git_version = check_output("git describe --tags --long --dirty --always".split()).decode('utf-8').strip()
    return format_version(version=git_version)


robbie_version = get_git_version()
#m ake a temporary version file to be installed then delete it
with open("robbie_version.sh", "a") as the_file:
    the_file.write(f"#!/bin/bash -l\necho {robbie_version}")

setup(
    name="Robbie",
    version=robbie_version,
    description="A batch processing work-flow for the detection of radio transients and variables",
    url="https://github.com/PaulHancock/Robbie",
    python_requires=">=3.6",
    scripts=[
        "robbie_version.sh",
        # python
        "scripts/auto_corr.py",
        "scripts/calc_var.py",
        "scripts/collect_transients.py",
        "scripts/filter_transients.py",
        "scripts/get_epoch.py",
        "scripts/get_lc_from_vot.py",
        "scripts/join_catalogues.py",
        "scripts/plot_variables.py",
        # nextflow
        "robbie.nf",
    ],
)
os.remove("robbie_version.sh")
