import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="proteo_estimator",
    version="0.0.5",
    author="Anna Calinawan",
    author_email="anna.calinawan@mssm.edu",
    description="A package from the NCI-CPTAC DREAM Proteogenomics Challenge",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Sage-Bionetworks/NCI-CPTAC-DREAM-Proteogenomics-Challenge/tree/master/Subchallenges_2_3/ProteoEstimator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
        "Operating System :: OS Independent",
    ],
    install_requires=[
       'docker>=4.0.1',
    ],
)