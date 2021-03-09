import os
import subprocess
import sys
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_call(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed")

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))

        cfg = 'Debug' if self.debug else 'Release'

        target_version = str(sys.version_info.major) + "." + str(
            sys.version_info.minor)

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DCMAKE_BUILD_TYPE=" + cfg,
            "-DPYBIND11_PYTHON_VERSION=" + target_version,
            "-DZIN_BUILD_APP=OFF",
            "-DZIN_BUILD_TESTS=OFF",
            "-DZIN_BUILD_PYTHON_BINDINGS=ON",
        ]

        cmake_list_dir = os.path.abspath(os.path.dirname(__file__))

        configure_command = ['cmake', cmake_list_dir] + cmake_args
        build_command = ['cmake', '--build', '.']

        subprocess.check_call(configure_command, cwd=self.build_temp)
        subprocess.check_call(build_command, cwd=self.build_temp)


setup(
    name="pysps",
    version="0.2",
    author="Yuki Koyama",
    author_email="yuki@koyama.xyz",
    description="Python bindings of sequential plane search",
    long_description="",
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("pysps")],
    url="",
    packages=find_packages(),
    install_requires=["numpy"],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.6",
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
