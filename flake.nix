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
    "aarch64-linux"
    "x86_64-darwin"
    "aarch64-darwin"
  ] (system:
    let
      pkgs = import nixpkgs {
        inherit system;
      };

      kmtricksBuildInputs = [
        pkgs.gcc13
        pkgs.zlib
        pkgs.cmake
      ];

      kmtricks = (with pkgs; stdenvNoCC.mkDerivation {
        pname = "kmtricks";
        version = "1.5.0";
        src = builtins.fetchGit {
          url = "https://github.com/tlemane/kmtricks";
          rev = "d1d37dd03beef497c143c8ed711980cd4c07b9b4";
          submodules = true;
        };

        nativeBuildInputs = kmtricksBuildInputs;

        configurePhase = ''
          cmake -S . -B build -DWITH_MODULES=ON
        '';

        buildPhase = ''
          cmake --build ./build --parallel 4
        '';

        installPhase = ''
          mkdir -p $out/bin
          cp ./bin/kmtricks $out/bin
        '';
      });

    in rec {
      defaultApp = flake-utils.lib.mkApp {
          drv = defaultPackage;
      };
      defaultPackage = kmtricks;
      devShell = pkgs.mkShellNoCC {
        name = "kmtricks_dev";
        buildInputs = kmtricksBuildInputs;
      };
      app.default = {
        type = "app";
        program = "${self.defaultPackage.${system}}/bin/kmtricks";
      };

    }
  );
}
