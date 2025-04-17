from conans import ConanFile
from conan.tools.cmake import CMakeDeps, CMake, CMakeToolchain
from conans.tools import save, load
import os
import pathlib
import subprocess
from rules_support import PluginBranchInfo


class CrossSpeciesComparisonGeneDetectPluginConan(ConanFile):
    """Class to package CrossSpeciesComparisonGeneDetectPlugin using conan

    Packages both RELEASE and RELWITHDEBINFO.
    Uses rules_support (github.com/ManiVaultStudio/rulessupport) to derive
    versioninfo based on the branch naming convention
    as described in https://github.com/ManiVaultStudio/core/wiki/Branch-naming-rules
    """

    name = "CrossSpeciesComparisonGeneDetectPlugin"
    description = (
        "CrossSpeciesComparisonGeneDetectPlugin."
    )
    topics = ("hdps", "plugin", "CrossSpeciesComparisonGeneDetectPlugin", "various")
    url = "https://github.com/ManiVaultStudio/CrossSpeciesComparisonGeneDetectPlugin"
    author = "B. van Lew b.van_lew@lumc.nl"  # conan recipe author
    license = "MIT"

    short_paths = True
    generators = "CMakeDeps"

    # Options may need to change depending on the packaged library
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
        path = load(
            pathlib.Path(pathlib.Path(__file__).parent.resolve(), "__gitpath.txt")
        )
        print(f"git info from {path}")
        return path

    def export(self):
        print("In export")
        # save the original source path to the directory used to build the package
        save(
            pathlib.Path(self.export_folder, "__gitpath.txt"),
            str(pathlib.Path(__file__).parent.resolve()),
        )

    def set_version(self):
        # Assign a version from the branch name
        branch_info = PluginBranchInfo(self.recipe_folder)
        self.version = branch_info.version

    def requirements(self):
        branch_info = PluginBranchInfo(self.__get_git_path())
        print(f"Core requirement {branch_info.core_requirement}")
        self.requires(branch_info.core_requirement)
        self.requires("CrossSpeciesComparisonTreeData/cytosploreviewer@lkeb/stable")

    # Remove runtime and use always default (MD/MDd)
    def configure(self):
        pass

    def system_requirements(self):
        #  May be needed for macOS or Linux
        pass

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def generate(self):
        generator = None
        if self.settings.os == "Macos":
            generator = "Xcode"
        if self.settings.os == "Linux":
            generator = "Ninja Multi-Config"

        tc = CMakeToolchain(self, generator=generator)

        # Use the Qt provided .cmake files
        qt_path = pathlib.Path(self.deps_cpp_info["qt"].rootpath)
        qt_cfg = list(qt_path.glob("**/Qt6Config.cmake"))[0]
        qt_dir = qt_cfg.parents[0].as_posix()

        tc.variables["Qt6_DIR"] = qt_dir

        # Use the ManiVault .cmake file to find ManiVault with find_package
        mv_core_root = self.deps_cpp_info["hdps-core"].rootpath
        manivault_dir = pathlib.Path(mv_core_root, "cmake", "mv").as_posix()
        print("ManiVault_DIR: ", manivault_dir)
        tc.variables["ManiVault_DIR"] = manivault_dir

        MV_CSCTD_PATH = pathlib.Path(self.deps_cpp_info["CrossSpeciesComparisonTreeData"].rootpath).as_posix()
        print(f"MV_CSCTD_INSTALL_DIR: {MV_CSCTD_PATH}")
        tc.variables["MV_INSTALL_DIR"] = self.install_dir
        tc.variables["MV_CSCTD_INSTALL_DIR"] = MV_CSCTD_PATH

        # Set some build options
        tc.variables["MV_UNITY_BUILD"] = "ON"

        # Use vcpkg-installed dependencies if there are any
        if os.environ.get("VCPKG_ROOT", None):
            vcpkg_dir = pathlib.Path(os.environ["VCPKG_ROOT"])
            vcpkg_exe = vcpkg_dir / "vcpkg.exe" if self.settings.os == "Windows" else vcpkg_dir / "vcpkg" 
            vcpkg_tc  = vcpkg_dir / "scripts" / "buildsystems" / "vcpkg.cmake"

            vcpkg_triplet = "x64-windows"
            if self.settings.os == "Macos":
                vcpkg_triplet = "x64-osx"
            if self.settings.os == "Linux":
                vcpkg_triplet = "x64-linux"

            print("vcpkg_dir: ", vcpkg_dir)
            print("vcpkg_exe: ", vcpkg_exe)
            print("vcpkg_tc: ", vcpkg_tc)
            print("vcpkg_triplet: ", vcpkg_triplet)

            tc.variables["VCPKG_LIBRARY_LINKAGE"]   = "dynamic"
            tc.variables["VCPKG_TARGET_TRIPLET"]    = vcpkg_triplet
            tc.variables["VCPKG_HOST_TRIPLET"]      = vcpkg_triplet
            tc.variables["VCPKG_ROOT"]              = vcpkg_dir.as_posix()

            tc.variables["CMAKE_PROJECT_INCLUDE"]   = vcpkg_tc.as_posix()

        tc.generate()

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure(build_script_folder="hdps/CrossSpeciesComparisonGeneDetectPlugin")
        cmake.verbose = True
        return cmake

    def build(self):
        print("Build OS is: ", self.settings.os)

        cmake = self._configure_cmake()
        cmake.build(build_type="RelWithDebInfo")
        cmake.build(build_type="Release")

    def package(self):
        package_dir = pathlib.Path(self.build_folder, "package")
        relWithDebInfo_dir = package_dir / "RelWithDebInfo"
        release_dir = package_dir / "Release"
        print("Packaging install dir: ", package_dir)
        subprocess.run(
            [
                "cmake",
                "--install",
                self.build_folder,
                "--config",
                "RelWithDebInfo",
                "--prefix",
                relWithDebInfo_dir,
            ]
        )
        subprocess.run(
            [
                "cmake",
                "--install",
                self.build_folder,
                "--config",
                "Release",
                "--prefix",
                release_dir,
            ]
        )
        self.copy(pattern="*", src=package_dir)

    def package_info(self):
        self.cpp_info.relwithdebinfo.libdirs = ["RelWithDebInfo/lib"]
        self.cpp_info.relwithdebinfo.bindirs = ["RelWithDebInfo/Plugins", "RelWithDebInfo"]
        self.cpp_info.relwithdebinfo.includedirs = ["RelWithDebInfo/include", "RelWithDebInfo"]
        self.cpp_info.release.libdirs = ["Release/lib"]
        self.cpp_info.release.bindirs = ["Release/Plugins", "Release"]
        self.cpp_info.release.includedirs = ["Release/include", "Release"]
