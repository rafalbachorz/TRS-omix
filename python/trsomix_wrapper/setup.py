from setuptools import setup

setup(
    name="trsomix_wrapper",
    version="0.2",
    #package_dir={"": "src"},
    #include_package_data=True,
    py_modules=["trsomix_wrapper"],
    setup_requires=["cffi>=1.0.0"],
    install_requires=["cffi>=1.0.0"],
    cffi_modules=["trsomix_wrapper.py:ffibuilder"],
)