from conans import ConanFile
from conan.tools.cmake import CMakeDeps, CMake, CMakeToolchain
from conans.tools import save, load
import os
import pathlib
import subprocess
from rules_support import PluginBranchInfo

class CrossSpeciesComparisonGeneDetectPluginConan(ConanFile):
    name = "CrossSpeciesComparisonGeneDetectPlugin"
    description = "A plugin for data heat-maps in ManiVault."
    topics = ("hdps", "ManiVault", "plugin", "CrossSpeciesComparisonGeneDetectPlugin", "data visualization")
    url = "https://github.com/ManiVaultStudio/CrossSpeciesComparisonGeneDetectPlugin"
    author = "B. van Lew b.van_lew@lumc.nl"
    license = "MIT"

    short_paths = True
    generators = "CMakeDeps"

    settings = {"os": None, "build_type": None, "compiler": None, "arch": None}
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": True, "fPIC": True}

    scm = {
        "type": "git",
        "subfolder": "hdps/CrossSpeciesComparisonGeneDetectPlugin",
        "url": "auto",
        "revision": "auto",
    }

    def __get_git_path(self):
        path = load(pathlib.Path(pathlib.Path(__file__).parent.resolve(), "__gitpath.txt"))
        print(f"git info from {path}")
        return path

    def export(self):
        print("In export")
        save(
            pathlib.Path(self.export_folder, "__gitpath.txt"),
            str(pathlib.Path(__file__).parent.resolve()),
        )

    def set_version(self):
        branch_info = PluginBranchInfo(self.recipe_folder)
        self.version = branch_info.version

    def requirements(self):
        branch_info = PluginBranchInfo(self.__get_git_path())
        print(f"Core requirement {branch_info.core_requirement}")
        self.requires(branch_info.core_requirement)
        # Add TreeData dependency here instead of at class level
        self.requires("CrossSpeciesComparisonTreeData/cytosploreviewer@lkeb/stable")

    def configure(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def generate(self):
        generator = None
        if self.settings.os == "Macos":
            generator = "Xcode"
        if self.settings.os == "Linux":
            generator = "Ninja Multi-Config"

        tc = CMakeToolchain(self, generator=generator)
        tc.variables["CMAKE_CXX_STANDARD_REQUIRED"] = "ON"

        # Qt6 configuration
        qt_path = pathlib.Path(self.deps_cpp_info["qt"].rootpath)
        qt_cfg = list(qt_path.glob("**/Qt6Config.cmake"))[0]
        tc.variables["Qt6_DIR"] = qt_cfg.parents[0].as_posix()

        # TreeData plugin configuration
        csctd_root = pathlib.Path(self.deps_cpp_info["CrossSpeciesComparisonTreeData"].rootpath)
        tc.variables["MV_CSCTD_INSTALL_DIR"] = csctd_root.as_posix()
        print(f"Set MV_CSCTD_INSTALL_DIR to: {csctd_root.as_posix()}")

        # ManiVault configuration
        mv_core_root = self.deps_cpp_info["hdps-core"].rootpath
        manivault_dir = pathlib.Path(mv_core_root, "cmake", "mv").as_posix()
        tc.variables["ManiVault_DIR"] = manivault_dir

        tc.generate()

        # Generate deps
        deps = CMakeDeps(self)
        deps.generate()

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(build_script_folder="hdps/CrossSpeciesComparisonGeneDetectPlugin")
        cmake.verbose = True
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build(build_type="RelWithDebInfo")
        cmake.build(build_type="Release")

    def package(self):
        cmake = CMake(self)
        cmake.install(build_type="RelWithDebInfo")
        cmake.install(build_type="Release")

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "both")
        self.cpp_info.set_property("cmake_file_name", "CrossSpeciesComparisonGeneDetectPlugin")
        
        # For RelWithDebInfo configuration
        self.cpp_info.components["_relwithdebinfo"].libdirs = ["RelWithDebInfo/Plugins"]
        self.cpp_info.components["_relwithdebinfo"].bindirs = ["RelWithDebInfo/Plugins"]
        self.cpp_info.components["_relwithdebinfo"].includedirs = ["RelWithDebInfo/include"]
        
        # For Release configuration
        self.cpp_info.components["_release"].libdirs = ["Release/Plugins"]
        self.cpp_info.components["_release"].bindirs = ["Release/Plugins"]
        self.cpp_info.components["_release"].includedirs = ["Release/include"]
