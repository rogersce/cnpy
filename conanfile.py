import os

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy, rmdir

required_conan_version = ">=2.0"


class CnpyConan(ConanFile):
    name = "cnpy"
    version = "2.0.0"
    license = "MIT"
    url = "https://github.com/lelapinagile/cnpy"
    homepage = "https://github.com/lelapinagile/cnpy"
    description = "Read and write NumPy .npy and .npz files from C++20"
    topics = ("numpy", "npy", "npz", "serialization", "scientific-computing")
    package_type = "library"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    exports_sources = (
        "CMakeLists.txt", "cmake/*", "include/*", "src/*", "tests/*", "examples/*", "bench/*", "LICENSE",
    )

    @property
    def _run_tests(self):
        return not self.conf.get("tools.build:skip_test", default=False, check_type=bool)

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        # zlib is an implementation detail: cnpy's public headers do not include zlib.h
        self.requires("zlib/[>=1.2.11 <2]")

    def build_requirements(self):
        if self._run_tests:
            self.test_requires("doctest/2.4.11")

    def validate(self):
        check_min_cppstd(self, 20)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CNPY_BUILD_TESTS"] = self._run_tests
        tc.cache_variables["CNPY_BUILD_EXAMPLES"] = False
        tc.cache_variables["CNPY_INSTALL"] = True
        tc.generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        if self._run_tests:
            cmake.test()

    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        CMake(self).install()
        rmdir(self, os.path.join(self.package_folder, "lib", "cmake"))  # CMakeDeps generates its own config

    def package_info(self):
        self.cpp_info.libs = ["cnpy"]
        self.cpp_info.set_property("cmake_file_name", "cnpy")
        self.cpp_info.set_property("cmake_target_name", "cnpy::cnpy")
        self.cpp_info.requires = ["zlib::zlib"]
        if not self.options.shared:
            self.cpp_info.defines = ["CNPY_STATIC_DEFINE"]
