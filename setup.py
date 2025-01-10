from setuptools import setup, find_packages
from pathlib import Path

long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name="burstlink",                  
    version="1.0.5",   
    description="A user-friendly package for analyzing gene interactions and transcriptional bursting.",                
    packages=find_packages(),            
    python_requires=">=3.8",   
    install_requires=[
        "statsmodels>=0.12.0",
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "pyarrow",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",  
)