

- set ONEAPI_ROOT in .bashrc
- run `source $ONEAPI_ROOT/setvars.sh` in .bashrc
- install make and cmake
- sudo apt install binutils
- sudo apt install g++ (and gcc)

- set up .vscode/c_cpp_properties.json (see file in example directory)

- C/C++ runner extension in VSCode

- add /snap/cmake/current/bin/cmake to cmake.cmakePath in settings
    - cmake on /snap/bin does not work with CMake tools in VSCode because it is a symlink
- ctrl-Shift P: Edit User-Local CMake kits
    - Add item:
    ```
    {
        "name": "Intel oneAPI",
        "compilers": {
        "C": "/s/intel/oneapi/compiler/2023.2.0/linux/bin/icx",
        "CXX": "/s/intel/oneapi/compiler/2023.2.0/linux/bin/icpx"
        },
        "isTrusted": true
    }
    ```
    - optional: delete the default gcc item
- ctrl-Shift P: Select CMake Variant
    - Can select Debug, Release, etc. from here
- add
    ```
    "cmake.debugConfig": {
        "MIMode": "gdb",
        "miDebuggerPath": "/s/intel/oneapi/debugger/2023.2.0/gdb/intel64/bin/gdb-oneapi"
    }
    ```
    to settings.json (can be accessed by finding the Debugger path setting in CMake tools extension settings)