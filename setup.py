import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pythonDMDapp",
    version="0.0.2",
    author="Daniel S Grün",
    author_email="danielgruns@gmail.com",
    description="An easy Vialux DMD package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/danielsgrun/pythonDMDapp",
    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["pythonDMDapp"],
    python_requires=">=3.6",
)