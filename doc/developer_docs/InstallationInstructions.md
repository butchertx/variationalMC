# Using variationalMC as a library

TBD

# Installing the full development environment

## Cloning the repository

Click the green "Code" button at the top right of the homepage for the repository.  Make sure the HTTPS tab is selected, copy the https address, then on your system navigate to the directory where you want to keep the code.  Type the command (in Git Bash for Windows, or terminal for Linux/Unix/Mac):

``` bash
git clone https://github.com/butchertx/variationalMC.git
```

This will copy the variationalMC directory into your current working directory.

## Install Intel oneAPI

Get the Intel oneAPI toolkit here: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html. Download the correct version for the platform you're working on. Development on this project is primarily supported for Linux/WSL, but the setup is intended to be flexible enough for cross-platform work. These instructions were written for the 2025.1.0 version of the Intel oneAPI Base Toolkit.

When you have downloaded the installer, follow the instructions on the Intel webpage (linked above) to install and configure your system. Once it is installed, the main action for configuration is adding the oneAPI configuration script to your .bashrc file. This will configure the environment variables that are needed for compiling using the Intel compilers. Add the following line at the end of the file `~/.bashrc`:

``` bash
source /opt/intel/oneapi/2025.1/oneapi-vars.sh
```
If using a different version of the oneAPI toolkit, be sure to replace the "2025.1" with whichever version is installed. You can view the installed versions with the command

``` bash
ls /opt/intel/oneapi
```

## Set up CMake

CMake is installed as part of the oneAPI configuration, so it should already be available on your system. Check this by running `which cmake`, and if it returns a path then cmake is already available. If it is not, run the following (Ubuntu):

``` bash
sudo apt update
sudo apt -y install cmake pkg-config build-essential
```

## Pull in Gtest for testing

- download the `googletest` repository from git: 
  - recommended method is navigate to https://github.com/google/googletest, click "Code" in the top right, and choose "Download Zip"
- On WSL:Ubuntu, I run the following commands to pull in the .zip from the default download location and unzip it:

```
cd ~/projects/variationalMC
mkdir -p thirdparty
cd thirdparty
cp /mnt/c/Users/Matthew/Downloads/googletest-main.zip .
unzip googletest-main.zip -d .
```
- The result should be the following directory structure: `thirdparty/googletest-main/...`, where `...` are the files in the root directory of the `googletest` repository.