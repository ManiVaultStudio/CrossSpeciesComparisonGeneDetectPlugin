from conan import ConanFile
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain
from conan.tools.files import load, save
from conan.tools.scm import Git
import os
import pathlib
import subprocess
from rules_support import PluginBranchInfo

class CrossSpeciesComparisonGeneDetectPluginConan(ConanFile):
    name = "CrossSpeciesComparisonGeneDetectPlugin"
    description = "A plugin for cross-species gene comparison in ManiVault"
    topics = ("hdps", "manivault", "plugin", "bioinformatics")
    url = "https://github.com/ManiVaultStudio/CrossSpeciesComparisonGeneDetectPlugin"
    author = "B. van Lew <b.van_lew@lumc.nl>"
    license = "MIT"

    # Build configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }
    default_options = {
        "shared": True,
        "fPIC": True,
    }

    # SCM configuration for in-source development
    scm = {
        "type": "git",
        "subfolder": ".",
        "url": "auto",
        "revision": "auto",
    }

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def set_version(self):
        branch_info = PluginBranchInfo(self.recipe_folder)
        self.version = branch_info.version

    def requirements(self):
        # Core dependency
        branch_info = PluginBranchInfo(self.recipe_folder)
        self.requires(branch_info.core_requirement)
        
        # Qt dependency
        self.requires("qt/6.8.2@lkeb/stable")
        
        # TreeData plugin dependency - using version range for flexibility
        self.requires("CrossSpeciesComparisonTreeData/cytosploreviewer@lkeb/stable")

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.22]")
        if self.settings.os == "Windows":
            self.tool_requires("ninja/[>=1.11.1]")

    def layout(self):
        self.folders.source = "."
        self.folders.build = os.path.join("build", str(self.settings.build_type))
        self.folders.generators = os.path.join(self.folders.build, "generators")

    def generate(self):
        # Generate CMake toolchain
        tc = CMakeToolchain(self)
        
        # Set CMake variables
        tc.variables["CMAKE_CXX_STANDARD"] = "20"
        tc.variables["CMAKE_CXX_STANDARD_REQUIRED"] = "ON"
        
        # Handle TreeData plugin path
        tree_data_path = pathlib.Path(self.dependencies["CrossSpeciesComparisonTreeData"].package_folder)
        tc.variables["MV_CSCTD_INSTALL_DIR"] = tree_data_path.as_posix()
        
        # Qt configuration
        qt_path = pathlib.Path(self.dependencies["qt"].package_folder)
        qt_cmake_path = next(qt_path.glob("**/Qt6Config.cmake")).parent
        tc.variables["Qt6_DIR"] = qt_cmake_path.as_posix()
        
        # ManiVault configuration
        mv_path = pathlib.Path(self.dependencies["hdps-core"].package_folder)
        tc.variables["ManiVault_DIR"] = mv_path.joinpath("cmake", "mv").as_posix()
        
        tc.generate()
        
        # Generate dependencies
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "CrossSpeciesComparisonGeneDetectPlugin")
        self.cpp_info.set_property("cmake_target_name", "ManiVault::CrossSpeciesComparisonGeneDetectPlugin")
        
        # Plugin binaries
        if self.settings.os == "Windows":
            self.cpp_info.bindirs = ["Plugins"]
            self.cpp_info.libdirs = ["Plugins"]
        else:
            self.cpp_info.libdirs = ["lib"]
        
        # Include directories
        self.cpp_info.includedirs = ["include"]
        
        # Runtime dependencies
        self.cpp_info.requires = [
            "hdps-core::hdps-core",
            "qt::qt",
            "CrossSpeciesComparisonTreeData::CrossSpeciesComparisonTreeData"
        ]
