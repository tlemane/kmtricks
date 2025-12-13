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

      kmtricks = (with pkgs; stdenv.mkDerivation {
        pname = "kmtricks";
        version = "1.5.0";
        src = builtins.fetchGit {
          url = "https://github.com/tlemane/kmtricks";
          rev = "561554f34e08a448e4f4a85ff3b0188ae938c9ae";
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
