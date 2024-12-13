from setuptools import setup, find_packages
from pathlib import Path

long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name="burstlink",                  
    version="0.3.6",                   
    packages=find_packages(),            
    python_requires=">=3.9",            
    long_description=long_description,
    long_description_content_type="text/markdown",  
)
