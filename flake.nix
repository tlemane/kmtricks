{
  inputs = {
    nixpkgs = {
      url = "github:nixos/nixpkgs/nixos-unstable";
    };
    flake-utils = {
      url = "github:numtide/flake-utils";
    };
  };
  outputs = { self, nixpkgs, flake-utils, ... }: flake-utils.lib.eachSystem [
    "x86_64-linux"
  ] (system:
    let
      pkgs = import nixpkgs {
        inherit system;
      };

      kmtricksBuildInputs = [
        pkgs.gcc9
        pkgs.cmake
        pkgs.zlib
        pkgs.gbenchmark
        pkgs.bzip2
        pkgs.gtest
        pkgs.zstd
        pkgs.virtualenv
      ];
    in rec {
      devShell = pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; } {
        name = "kmtricks_dev_env";
        buildInputs = kmtricksBuildInputs;
      };
    }
  );
}
