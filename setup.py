from setuptools import setup, find_packages, Extension

module = Extension(
            "spkmeansmodule",
            sources=['spkmeans.c', "spkmeansmodule.c"])


setup(
    name='spkmeansmodule',
    version='0.1.0',
    author='Yaakov Goldsmith and Lilach Daniel final project',
    description="spkmeans",
    install_requires=["invoke"],
    packages=find_packages(),
    license='GPL-2',
    classifiers=['Development Status :: 3 - Alpha'],
    ext_modules=[module]

)
