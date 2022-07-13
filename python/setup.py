from setuptools import setup

setup(
    name="trsomix",
    version="0.1",
    py_modules=["trsomix"],
    setup_requires=["cffi"],
    install_requires=["cffi"],
    cffi_modules=["build_trsomix:ffivuilder"]
)