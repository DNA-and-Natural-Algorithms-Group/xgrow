{ pkgs ? import <nixpkgs> {}, pythonPackages ? pkgs.python3Packages }:

pythonPackages.buildPythonPackage {
  name = "xgrow";
  src = ./.;
  propagatedBuildInputs = [ pkgs.xorg.libX11
                            pythonPackages.numpy
                            pythonPackages.pandas
                            pythonPackages.matplotlib
                            pythonPackages.pyyaml ];
}
