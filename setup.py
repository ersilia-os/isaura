import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()

setuptools.setup(
    name="isaura",
    version="0.0.1",
    author="Ersilia Open Source Initiative",
    author_email="hello@ersilia.io",
    description="Isaura data lake for pre-computed Ersilia properties",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ersilia-os/isaura",
    project_urls={"GitBook": "https://ersilia.gitbook.io/ersilia/",},
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(exclude=["utilities"]),
    python_requires=">=3.7",
)
