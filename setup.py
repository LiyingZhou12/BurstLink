from setuptools import setup, find_packages
from pathlib import Path

long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name="burstlink",                  
    version="0.0.2",   
    description="A user-friendly package for analyzing gene interactions and transcriptional bursting.",                
    packages=find_packages(),            
    python_requires=">=3.8.18",            
    long_description=long_description,
    long_description_content_type="text/markdown",  
)
