

Installation guide
################################

Installing the full development environment
===========================================

Cloning the repository
-------------------------------------------

Click the green "Code" button at the top right of the homepage for the repository.  Make sure the HTTPS tab is selected, copy the https address, then on your system navigate to the directory where you want to keep the code.  Type the command (in Git Bash for Windows, or terminal for Linux/Unix/Mac):

.. code-block:: bash

    git clone https://github.com/butchertx/variationalMC.git


This will copy the variationalMC directory into your current working directory.

Install Intel oneAPI
-------------------------------------------

Get the Intel oneAPI toolkit [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). Download the correct version for the platform you're working on. Development on this project is primarily supported for Linux/WSL, but the setup is intended to be flexible enough for cross-platform work. These instructions were written for the 2025.1.0 version of the Intel oneAPI Base Toolkit built on Ubuntu-22.04 within WSL.

When you have downloaded the installer, follow the instructions on the Intel webpage (linked above) to install and configure your system. Once it is installed, the main action for configuration is adding the oneAPI configuration script to your .bashrc file. This will configure the environment variables that are needed for compiling using the Intel compilers. Add the following line at the end of the file ``~/.bashrc``:

.. code-block:: bash
    
    source /opt/intel/oneapi/2025.1/oneapi-vars.sh

If using a different version of the oneAPI toolkit, be sure to replace the "2025.1" with whichever version is installed. You can view the installed versions with the command

.. code-block:: bash
    
    ls /opt/intel/oneapi

Install the build dependencies via conda
-------------------------------------------

The build system uses ``cmake``, and the documentation is built using sphinx and doxygen. To install all of the build dependencies together, it is easiest to use conda. If you have conda installed, the recommended way to do this is to create a new conda environment. This conda environment is specified in the ``doc/conda.yml`` file (from the root directory of the repository). To create the conda environment, run the following command from the root of the repository in the terminal:

.. code-block:: bash
    
    conda env create -f doc/conda.yml

Then, activate the environment with:

.. code-block:: bash
    
    conda activate vmc_build

If this doesn't work, you can also try building an environment manually using ``doc/requirements.txt``.

Pull in Gtest for testing
-------------------------------------------

* download the ``googletest`` repository from git:
  * recommended method is navigate to the [googletest Github repository](https://github.com/google/googletest), click "Code" in the top right, and choose "Download Zip"
* On WSL:Ubuntu, I run the following commands to pull in the .zip from the default download location and unzip it:

.. code-block:: bash
    
    cd ~/projects/variationalMC
    mkdir -p thirdparty
    cd thirdparty
    cp /mnt/c/Users/Matthew/Downloads/googletest-main.zip .
    unzip googletest-main.zip -d .

* The result should be the following directory structure: ``thirdparty/googletest-main/...``, where ``...`` are the files in the root directory of the ``googletest`` repository.

Optional: Set up VS Code as your IDE
-------------------------------------------

This is the recommended IDE setup because it comes with extension support for CMake, testing, and C++ language features. It also has a RemoteSSH extension that will allow you to work through the same interface on a remote host (compute cluster). With the ``variationalMC`` directory open in VS Code, do the following:

Install the following extensions:

* C/C++
* CMake Tools

Configuration
+++++++++++++++++++++++++++++++++++++++++++

#. Configure the CMake extension settings. You can get the ``cmake`` path by typing ``which cmake`` in the terminal. If ``.vscode/settings.json`` does not exist, create it and populate it with the following:

    .. code-block:: json
        
        {
            "cmake.cmakePath": "/usr/bin/cmake",
            "cmake.configureOnOpen": true,
            "cmake.buildDirectory": "${workspaceFolder}/build",
            "cmake.debugConfig": {
                "MIMode": "gdb",
                "miDebuggerPath": " /opt/intel/oneapi/debugger/2025.1/bin/gdb-oneapi"
            }
        }

    If the ``settings.json`` file already exists, just add the above options alongside any other options already in there.

#. Press Ctrl+Shift+P and Select "Edit User-Local CMake kits" (start typing it and VSCode will search for it). Add the following configuration to the resulting .json file:

    .. code-block:: json
        
        {
            "name": "Intel oneAPI",
            "compilers": {
                "C": "/opt/intel/oneapi/compiler/2025.1/bin/icx",
                "CXX": "/opt/intel/oneapi/compiler/2025.1/bin/icpx"
            },
            "isTrusted": true
        }

    Optional: delete the default configuration(s) that was already in this file.


(Optional) Set up SSH for remote development
-------------------------------------------
