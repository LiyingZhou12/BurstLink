from setuptools import setup, find_packages
from pathlib import Path

long_description = Path("README.md").read_text(encoding="utf-8")

def parse_requirements(filename):
    with open(filename, "r") as file:
        return [line.strip() for line in file if line and not line.startswith("#")]

setup(
    name="burstlink",                  
    version="1.0.0",   
    description="A user-friendly package for analyzing gene interactions and transcriptional bursting.",                
    packages=find_packages(),            
    python_requires=">=3.8",   
    install_requires=parse_requirements("requirements.txt"),
    long_description=long_description,
    long_description_content_type="text/markdown",  
)
